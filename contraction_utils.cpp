#include "contraction_utils.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <regex>
#include <sstream>
#include <vector>
#include <unordered_set>

namespace qflex {

// ContractionData methods

ContractionData ContractionData::Initialize(
    const std::list<ContractionOperation>& ordering,
    std::vector<std::vector<MKLTensor>>* tensor_grid,
    std::vector<std::complex<double>>* amplitudes) {
  ContractionData data;
  data.tensor_grid_ = tensor_grid;
  data.amplitudes_ = amplitudes;

  // Maximum bond dimension in the tensor grid.
  int bond_dim = 0;

  // Indices of tensor_grid elements during the contraction.
  std::vector<std::vector<std::vector<std::string>>> grid_indices(
      tensor_grid->size());
  for (int i = 0; i < tensor_grid->size(); ++i) {
    grid_indices[i] =
        std::vector<std::vector<std::string>>((*tensor_grid)[i].size());
    for (int j = 0; j < (*tensor_grid)[i].size(); ++j) {
      grid_indices[i][j] = (*tensor_grid)[i][j].get_indices();
      for (const int dim : (*tensor_grid)[i][j].get_dimensions()) {
        bond_dim = std::max(bond_dim, dim);
      }
    }
  }

  // Indices of each patch.
  std::unordered_map<std::string, std::vector<std::string>> patch_indices;
  // Rank of cut-copy tensors.
  std::unordered_map<std::string, int> cut_copy_rank;

  for (const auto& op : ordering) {
    switch (op.op_type) {
      case ContractionOperation::EXPAND: {
        const auto& indices_a = patch_indices[op.expand.id];
        const auto& indices_b =
            grid_indices[op.expand.tensor[0]][op.expand.tensor[1]];
        patch_indices[op.expand.id] =
            _vector_subtraction(_vector_union(indices_a, indices_b),
                                _vector_intersection(indices_a, indices_b));
        if (data.patch_rank_[op.expand.id] <
            patch_indices[op.expand.id].size()) {
          data.patch_rank_[op.expand.id] = patch_indices[op.expand.id].size();
        }
        break;
      }
      case ContractionOperation::CUT: {
        const std::string cut_index = index_name(op.cut.tensors);
        int side = 0;
        for (auto& tensor : op.cut.tensors) {
          auto& indices = grid_indices[tensor[0]][tensor[1]];
          const std::string copy_name = cut_copy_name(op.cut.tensors, side);
          cut_copy_rank[copy_name] = indices.size();
          ++side;
          indices = _vector_subtraction(indices, {cut_index});
        }
        break;
      }
      case ContractionOperation::MERGE: {
        const auto& indices_a = patch_indices[op.merge.target_id];
        const auto& indices_b = patch_indices[op.merge.source_id];
        patch_indices[op.merge.target_id] =
            _vector_subtraction(_vector_union(indices_a, indices_b),
                                _vector_intersection(indices_a, indices_b));
        if (data.patch_rank_[op.merge.target_id] <
            patch_indices[op.merge.target_id].size()) {
          data.patch_rank_[op.merge.target_id] =
              patch_indices[op.merge.target_id].size();
        }
        break;
      }
    }
  }

  long allocated_space = 0;

  // Max rank/size of all patches.
  int max_rank = 0;
  for (const auto& patch_rank_pair : data.patch_rank_) {
    max_rank = std::max(patch_rank_pair.second, max_rank);
  }
  int max_size = (int)pow(bond_dim, max_rank);

  // General-purpose scratch space (primarily used for tensor reordering).
  data.scratch_.push_back(MKLTensor({""}, {max_size}));
  data.scratch_map_[kGeneralSpace] = 0;
  allocated_space += max_size;

  // "Swap tensor" space, used to store operation results.
  for (int rank = 1; rank <= max_rank; ++rank) {
    const int size = (int)pow(bond_dim, rank);
    data.scratch_.push_back(MKLTensor({""}, {size}));
    data.scratch_map_[result_space(rank)] = rank;
    allocated_space += size;
  }

  int patch_pos = data.scratch_map_.size();
  for (const auto& patch_rank_pair : data.patch_rank_) {
    const int size = (int)pow(bond_dim, patch_rank_pair.second);
    data.scratch_.push_back(MKLTensor({""}, {size}));
    data.scratch_map_[patch_rank_pair.first] = patch_pos++;
    allocated_space += size;
  }

  // TODO(martinop): minor optimizations possible: When consecutive cuts apply
  // to the same grid tensor, only one copy needs to be stored.
  int cut_copy_pos = data.scratch_map_.size();
  for (const auto& copy_rank_pair : cut_copy_rank) {
    const int size = (int)pow(bond_dim, copy_rank_pair.second);
    data.scratch_.push_back(MKLTensor({""}, {size}));
    data.scratch_map_[copy_rank_pair.first] = patch_pos++;
    allocated_space += size;
  }

  // Log how much space is allocated for this operation. This includes all
  // tensors allocated during contraction; total memory usage should not vary
  // significantly from this value.
  int scale = 0;
  double alloc_size = allocated_space * sizeof(std::complex<double>);
  while (alloc_size > (1 << 10)) {
    ++scale;
    alloc_size /= (1 << 10);
  }
  int old_precision = std::cout.precision(6);
  std::string suffix[] = {"B", "kB", "MB", "GB"};
  std::cout << alloc_size << suffix[scale] << " allocated." << std::endl;
  std::cout.precision(old_precision);
  return data;
}

void ContractionData::ContractGrid(
    std::list<ContractionOperation> ordering, int output_index,
    std::unordered_map<std::string, bool> active_patches) {
  MKLTensor* output;
  while (!ordering.empty()) {
    const ContractionOperation op = ordering.front();
    ordering.pop_front();
    switch (op.op_type) {
      case ContractionOperation::EXPAND: {
        // Multiply by the new tensor and store result in scratch.
        MKLTensor& next =
            (*tensor_grid_)[op.expand.tensor[0]][op.expand.tensor[1]];
        if (!active_patches[op.expand.id]) {
          // First tensor in patch.
          get_scratch(op.expand.id) = next;
          active_patches[op.expand.id] = true;
          continue;
        }
        MKLTensor& prev = get_scratch(op.expand.id);
        std::string result_id = result_space(patch_rank_[op.expand.id]);
        MKLTensor& result = get_scratch(result_id);
        multiply(prev, next, result, get_scratch(kGeneralSpace).data());
        if (ordering.empty()) {
          output = &result;
          continue;
        }
        int temp = scratch_map_[op.expand.id];
        scratch_map_[op.expand.id] = scratch_map_[result_id];
        scratch_map_[result_id] = temp;
        continue;
      }
      case ContractionOperation::CUT: {
        const std::string index = index_name(op.cut.tensors);
        MKLTensor& tensor_a =
            (*tensor_grid_)[op.cut.tensors[0][0]][op.cut.tensors[0][1]];
        MKLTensor& copy_a =
            get_scratch(cut_copy_name(op.cut.tensors, /*side = */ 0));
        copy_a = tensor_a;
        // List of values to evaluate on the cut.
        std::vector<int> values = op.cut.values;
        if (values.empty()) {
          int num_values = tensor_a.get_index_to_dimension().at(index);
          for (int i = 0; i < num_values; ++i) {
            values.push_back(i);
          }
        }
        if (op.cut.tensors.size() > 1) {
          // This is a normal cut; each value adds to the same amplitude.
          MKLTensor& tensor_b =
              (*tensor_grid_)[op.cut.tensors[1][0]][op.cut.tensors[1][1]];
          MKLTensor& copy_b =
              get_scratch(cut_copy_name(op.cut.tensors, /*side = */ 1));
          copy_b = tensor_b;
          for (int val : values) {
            copy_a.project(index, val, tensor_a);
            copy_b.project(index, val, tensor_b);
            ContractGrid(ordering, output_index, active_patches);
          }
          tensor_a = copy_a;
          tensor_b = copy_b;
        } else {
          // This is a terminal cut; each value adds to a different amplitude.
          output_index *= values.size();
          for (int val : values) {
            copy_a.project(index, val, tensor_a);
            ContractGrid(ordering, output_index, active_patches);
            output_index++;
          }
          tensor_a = copy_a;
        }
        return;  // Post-cut contraction is handled recursively.
      }
      case ContractionOperation::MERGE: {
        // Multiply two existing tensors and store result in scratch.
        MKLTensor& patch_1 = get_scratch(op.merge.source_id);
        MKLTensor& patch_2 = get_scratch(op.merge.target_id);
        if (!active_patches[op.merge.target_id]) {
          // Copy the old patch into the new space.
          patch_2 = patch_1;
          active_patches[op.merge.target_id] = true;
          continue;
        }
        std::string result_id = result_space(patch_rank_[op.merge.target_id]);
        MKLTensor& result = get_scratch(result_id);
        multiply(patch_1, patch_2, result, get_scratch(kGeneralSpace).data());
        if (ordering.empty()) {
          output = &result;
          continue;
        }
        int temp = scratch_map_[op.merge.target_id];
        scratch_map_[op.merge.target_id] = scratch_map_[result_id];
        scratch_map_[result_id] = temp;
        continue;
      }
    }
  }
  // Check for an output of size larger than 1
  if (output->size() != 1) {
    std::cout << "Contraction did not complete; final tensor is ";
    output->print();
    assert(false);
  }
  (*amplitudes_)[output_index] += (*output->data());
  return;
}

// External methods

bool ordering_data_to_contraction_ordering(
    std::istream* ordering_data, const int I, const int J,
    const std::optional<std::vector<std::vector<int>>>& off,
    std::list<ContractionOperation>* ordering) {
  static const std::regex cut_value_regex("\\([0-9,]*\\)");
  std::string line;
  std::string error_msg;
  std::string operation;
  // Read one line at a time from the ordering, skipping comments.
  while (getline(*ordering_data, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::stringstream ss(line);
    // The first element is the operation (expand, cut, or merge).
    ss >> operation;
    if (operation == "expand") {
      std::string patch;
      int index;
      ss >> patch;
      ss >> index;
      if (index < 0) {
        error_msg = "Index cannot be negative.";
        break;
      }
      if (index >= I * J) {
        error_msg = "Index must be within grid boundaries.";
        break;
      }
      std::vector<int> position = get_qubit_coords(index, J);
      if (find_grid_coord_in_list(off, position[0], position[1])) {
        error_msg = "Index must specify an active qubit.";
        break;
      }
      ordering->emplace_back(ExpandPatch(patch, position));
    } else if (operation == "cut") {
      std::vector<int> values;
      std::string values_str;
      std::smatch match;
      ss >> values_str;
      if (!std::regex_match(values_str, match, cut_value_regex)) {
        error_msg = "Cut values must be comma-separated ints, e.g. (0,1,3).";
        break;
      }
      values_str = values_str.substr(1, values_str.size() - 2);
      if (!values_str.empty()) {
        std::string val;
        auto start = 0;
        auto end = values_str.find(',');
        while (end != std::string::npos) {
          values.push_back(atoi(values_str.substr(start, end).c_str()));
          start = end + 1;
          end = values_str.find(',', start);
        }
        values.push_back(atoi(values_str.substr(start, end).c_str()));
      }
      int index_1, index_2;
      ss >> index_1;
      if (index_1 < 0) {
        error_msg = "Index 1 cannot be negative.";
        break;
      }
      if (index_1 >= I * J) {
        error_msg = "Index 1 must be within grid boundaries.";
        break;
      }
      std::vector<int> position_1 = get_qubit_coords(index_1, J);
      if (find_grid_coord_in_list(off, position_1[0], position_1[1])) {
        error_msg = "Index 1 must specify an active qubit.";
        break;
      }
      if (ss.eof()) {
        ordering->emplace_back(CutIndex({position_1}, values));
      } else {
        ss >> index_2;
        if (index_2 < 0) {
          error_msg = "Index 2 cannot be negative";
          break;
        }
        if (index_2 >= I * J) {
          error_msg = "Index 2 must be within grid boundaries.";
          break;
        }
        std::vector<int> position_2 = get_qubit_coords(index_2, J);
        if (find_grid_coord_in_list(off, position_2[0], position_2[1])) {
          error_msg = "Index 2 must specify an active qubit.";
          break;
        }
        ordering->emplace_back(CutIndex({position_1, position_2}, values));
      }
    } else if (operation == "merge") {
      std::string patch_1, patch_2;
      ss >> patch_1;
      ss >> patch_2;
      ordering->emplace_back(MergePatches(patch_1, patch_2));
    } else {
      error_msg = "Received an invalid operation in config.";
      break;
    }
  }
  if (!error_msg.empty()) {
    std::cout << "Parsing failed on line: \"" << line
              << "\" with error: " << error_msg << std::endl;
    return false;
  }
  return true;
}

std::string index_name(const std::vector<int>& p1, const std::vector<int>& p2) {
  char buffer[64];
  if (p1.size() == 2 && p2.size() == 2) {
    // Two-qubit contraction.
    int len = snprintf(buffer, sizeof(buffer), "(%d,%d),(%d,%d)", p1[0], p1[1],
                       p2[0], p2[1]);
    return std::string(buffer, len);
  }
  if (p1.size() == 3 && p2.size() == 3) {
    // Single-qubit contraction, or virtual index.
    int len = snprintf(buffer, sizeof(buffer), "(%d,%d,%d),(%d,%d,%d)", p1[0],
                       p1[1], p1[2], p2[0], p2[1], p2[2]);
    return std::string(buffer, len);
  }
  // Final qubit output value assignment.
  if (p1.size() == 2 && p2.empty()) {
    int len = snprintf(buffer, sizeof(buffer), "(%d,%d),(o)", p1[0], p1[1]);
    return std::string(buffer, len);
  }
  // use for loop to print out p1 and p2
  std::cout << "Failed to construct tensor name with the following vectors: v1 = {";
  for(int i = 0; i < p1.size(); ++i) {
    std::cout << p1.at(i);
    if (i != p1.size() - 1) {
      std::cout << ", ";
    }
  }
  std::cout << "}, v2 = {";
  for(int j = 0; j < p2.size(); ) {
    std::cout << p2.at(j);
    if (j != p2.size() - 1) {
      std::cout << ", ";
    }
  }
  std::cout << "}." << std::endl;
  assert(false);
  return "";
}

std::string index_name(const std::vector<std::vector<int>>& tensors) {
  if (tensors.size() == 2) {
    return index_name(tensors.at(0), tensors.at(1));
  }
  if (tensors.size() == 1) {
    return index_name(tensors.at(0), {});
  }
  std::cout << "Failed to construct tensor name with input tensors size: " << tensors.size() << std::endl;
  assert(false);
  return "";
}

std::vector<int> get_qubit_coords(int q, int J) {
  int i = q / J;
  return std::vector<int>({i, q - i * J});
}

bool find_grid_coord_in_list(
    const std::optional<std::vector<std::vector<int>>>& coord_list, const int i,
    const int j) {
  return coord_list.has_value() &&
         find(coord_list.value().begin(), coord_list.value().end(),
              std::vector<int>({i, j})) != coord_list.value().end();
}

bool IsOrderingValid(const std::list<ContractionOperation>& ordering) {
  struct PatchState {
    // This patch has expanded, but no cuts have happened since then.
    bool is_active = false;
    // This patch has expanded and at least one cut happened since then.
    bool is_used = false;
    // This patch has been merged into another patch.
    bool is_merged = false;
  };
  std::unordered_map<std::string, PatchState> patches;
  std::unordered_set<std::string> cut_indices;
  std::unordered_set<std::string> used_tensors;
  int error_space = 4000;                  // ~4kiB for error logs
  char error_msg[error_space + 100] = "";  // space for "too long" warning
  int next_error = 0;                      // position of next error log
  for (const auto& op : ordering) {
    if (next_error >= error_space) {
      // The error log is too full to report any more errors; return early.
      next_error +=
          snprintf(error_msg + next_error, sizeof(error_msg) - next_error,
                   "Too many errors to log!\n");
      break;
    }
    switch (op.op_type) {
      case ContractionOperation::EXPAND: {
        if (patches[op.expand.id].is_used) {
          next_error += snprintf(
              error_msg + next_error, error_space - next_error,
              "Tensor at (%d,%d) is added to non-empty patch %s after a cut.\n",
              op.expand.tensor[0], op.expand.tensor[1], op.expand.id.c_str());
        }
        if (patches[op.expand.id].is_merged) {
          next_error += snprintf(
              error_msg + next_error, error_space - next_error,
              "Tensor at (%d,%d) is added to previously-merged patch %s.\n",
              op.expand.tensor[0], op.expand.tensor[1], op.expand.id.c_str());
        }
        char tensor_name[20];
        snprintf(tensor_name, sizeof(tensor_name), "(%d,%d)",
                 op.expand.tensor[0], op.expand.tensor[1]);
        if (used_tensors.find(tensor_name) != used_tensors.end()) {
          next_error += snprintf(
              error_msg + next_error, error_space - next_error,
              "Tensor %s is contracted multiple times.\n", tensor_name);
        }
        used_tensors.insert(tensor_name);
        patches[op.expand.id].is_active = true;
        continue;
      }
      case ContractionOperation::CUT: {
        for (auto& patch_pair : patches) {
          if (patch_pair.second.is_active) {
            patch_pair.second.is_used = true;
          }
        }
        const std::string index = index_name(op.cut.tensors);
        if (cut_indices.find(index) != cut_indices.end()) {
          next_error +=
              snprintf(error_msg + next_error, error_space - next_error,
                       "Index %s is cut multiple times.\n", index.c_str());
        }
        cut_indices.insert(index);
        continue;
      }
      case ContractionOperation::MERGE: {
        if (patches[op.merge.source_id].is_merged) {
          next_error +=
              snprintf(error_msg + next_error, error_space - next_error,
                       "Patch %s is merged multiple times.\n",
                       op.merge.source_id.c_str());
        }
        if (patches[op.merge.target_id].is_used) {
          next_error += snprintf(
              error_msg + next_error, error_space - next_error,
              "Patch %s is merged into non-empty patch %s after a cut.\n",
              op.merge.source_id.c_str(), op.merge.target_id.c_str());
        }
        patches[op.merge.source_id].is_merged = true;
        patches[op.merge.target_id].is_active = true;
        continue;
      }
    }
  }
  if (next_error == 0) {
    return true;
  }
  std::cout << error_msg << std::endl;
  return false;
}

void ContractGrid(const std::list<ContractionOperation>& ordering,
                  std::vector<std::vector<MKLTensor>>* tensor_grid,
                  std::vector<std::complex<double>>* amplitudes) {
  if (amplitudes == nullptr) {
    std::cout << "Amplitude return vector must be non-null." << std::endl;
    assert(amplitudes != nullptr);
  }
  if (IsOrderingValid(ordering) == false) {
    std::cout << "Ordering must be valid." << std::endl;
    assert(IsOrderingValid(ordering));
  }

  // Populate ContractionData and perform grid contraction.
  ContractionData data =
      ContractionData::Initialize(ordering, tensor_grid, amplitudes);
  std::unordered_map<std::string, bool> active_patches;
  for (const auto& patch : data.scratch_list()) {
    active_patches[patch] = false;
  }
  data.ContractGrid(ordering, /*output_index = */ 0, active_patches);
}

}  // namespace qflex
