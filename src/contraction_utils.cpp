#include "contraction_utils.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <unordered_set>

#include "errors.h"

namespace qflex {

// ContractionData methods

ContractionData ContractionData::Initialize(
    const std::list<ContractionOperation>& ordering,
    std::vector<std::vector<Tensor>>* tensor_grid,
    std::vector<std::complex<double>>* amplitudes) {
  if (tensor_grid == nullptr) {
    throw ERROR_MSG("Tensor grid must be non-null.");
  }
  if (amplitudes == nullptr) {
    throw ERROR_MSG("Amplitude return vector must be non-null.");
  }
  ContractionData data;
  data.tensor_grid_ = tensor_grid;
  data.amplitudes_ = amplitudes;

  // Map of index name to dimenion for all indices of all tensors  in the
  // network.
  std::unordered_map<std::string, std::size_t> index_dimension;

  // Indices of tensor_grid elements during the contraction.
  std::vector<std::vector<std::vector<std::string>>> grid_indices(
      tensor_grid->size());
  for (std::size_t i = 0; i < tensor_grid->size(); ++i) {
    grid_indices[i] =
        std::vector<std::vector<std::string>>((*tensor_grid)[i].size());
    for (std::size_t j = 0; j < (*tensor_grid)[i].size(); ++j) {
      grid_indices[i][j] = (*tensor_grid)[i][j].get_indices();
      for (const std::string& index : (*tensor_grid)[i][j].get_indices()) {
        index_dimension[index] =
            (*tensor_grid)[i][j].get_index_to_dimension().at(index);
      }
    }
  }

  // Indices of each patch.
  std::unordered_map<std::string, std::vector<std::string>> patch_indices;
  // Size of each patch.
  std::unordered_map<std::string, std::size_t> patch_size;
  // Size of cut-copy tensors.
  std::unordered_map<std::string, std::size_t> cut_copy_size;

  for (const auto& op : ordering) {
    switch (op.op_type) {
      case ContractionOperation::EXPAND: {
        const auto& indices_a = patch_indices[op.expand.id];
        const auto& indices_b =
            grid_indices[op.expand.tensor[0]][op.expand.tensor[1]];
        patch_indices[op.expand.id] =
            _vector_subtraction(_vector_union(indices_a, indices_b),
                                _vector_intersection(indices_a, indices_b));
        std::size_t new_size = 1;
        for (const auto& index : patch_indices[op.expand.id]) {
          new_size *= index_dimension.at(index);
        }
        patch_size[op.expand.id] = new_size;
        if (data.patch_max_size_[op.expand.id] < patch_size[op.expand.id]) {
          data.patch_max_size_[op.expand.id] = patch_size[op.expand.id];
        }
        break;
      }
      case ContractionOperation::CUT: {
        try {
          const std::string cut_index = index_name(op.cut.tensors);
          int side = 0;
          for (auto& tensor : op.cut.tensors) {
            auto& indices = grid_indices[tensor[0]][tensor[1]];
            const std::string copy_name = cut_copy_name(op.cut.tensors, side);
            ++side;
            // Update indices so it keeps track of cuts so that tensor sizes and
            // ranks are not overestimated.
            indices = _vector_subtraction(indices, {cut_index});
            // TODO(benjaminvillalonga): reduce size of copy, since it loses
            // an index. This has implications in ContractGrid() below.
            cut_copy_size[copy_name] =
                (*tensor_grid)[tensor[0]][tensor[1]].size();
          }
        } catch (const std::string& err_msg) {
          throw ERROR_MSG("Failed during CUT. Error:\n\t[", err_msg, "]");
        }
        break;
      }
      case ContractionOperation::MERGE: {
        const auto& indices_a = patch_indices[op.merge.target_id];
        const auto& indices_b = patch_indices[op.merge.source_id];
        patch_indices[op.merge.target_id] =
            _vector_subtraction(_vector_union(indices_a, indices_b),
                                _vector_intersection(indices_a, indices_b));
        std::size_t new_size = 1;
        for (const auto& index : patch_indices[op.merge.target_id]) {
          new_size *= index_dimension.at(index);
        }
        patch_size[op.merge.target_id] = new_size;
        if (data.patch_max_size_[op.merge.target_id] <
            patch_size[op.merge.target_id]) {
          data.patch_max_size_[op.merge.target_id] =
              patch_size[op.merge.target_id];
        }
        break;
      }
    }
  }

  // Calculate the necessary space before actually allocating memory.

  // If the --memory flag is set properly, this should prevent OOM errors.
  std::size_t allocated_space = 0;

  // Max size of a tensor (or patches) along the contraction.
  std::size_t max_size = 0;
  // The set of max_sizes of all patches.
  std::unordered_set<std::size_t> unique_sizes;
  for (const auto& patch_max_size_pair : data.patch_max_size_) {
    max_size = std::max(patch_max_size_pair.second, max_size);
    unique_sizes.insert(patch_max_size_pair.second);
  }

  // General-purpose scratch space (primarily used for tensor reordering).
  allocated_space += max_size;

  // "Swap tensor" space, used to store operation results.
  for (const auto& swap_size : unique_sizes) {
    allocated_space += swap_size;
  }
  // Per-patch space, sized to match the maximum size of the patch.
  for (const auto& patch_max_size_pair : data.patch_max_size_) {
    allocated_space += patch_max_size_pair.second;
  }
  // Extra space for storing copies made during cut operations.
  for (const auto& cut_copy_size_pair : cut_copy_size) {
    allocated_space += cut_copy_size_pair.second;
  }

  // Prevent memory allocation from exceeding memory_limit.
  double alloc_size = allocated_space * sizeof(s_type);
  if (alloc_size > global::memory_limit) {
    throw ERROR_MSG("Required space (", utils::readable_memory_string(alloc_size),
                    ") exceeds memory limit (",
                    utils::readable_memory_string(global::memory_limit),
                    "). Cancelling simulation.");
  }

  // Log how much space will be allocated for this operation. This includes all
  // tensors allocated during contraction; total memory usage should not vary
  // significantly from this value.
  if (global::verbose > 0) {
    std::size_t old_precision = std::cerr.precision(6);
    std::cerr << "Allocating " << utils::readable_memory_string(alloc_size)
              << " for this simulation." << std::endl;
    std::cerr.precision(old_precision);
  }

  // Actually allocate the required space.

  // General-purpose scratch space (primarily used for tensor reordering).
  try {
    data.scratch_.push_back(Tensor({""}, {max_size}));
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call Tensor(). Error:\n\t[", err_msg, "]");
  }
  data.scratch_map_[kGeneralSpace] = 0;
  // "Swap tensor" space, used to store operation results.
  {
    std::size_t patch_pos = 1;
    for (const auto& swap_size : unique_sizes) {
      try {
        data.scratch_.push_back(Tensor({""}, {swap_size}));
      } catch (const std::string& err_msg) {
        throw ERROR_MSG("Failed to call Tensor(). Error:\n\t[", err_msg, "]");
      }
      data.scratch_map_[result_space(swap_size)] = patch_pos;
      patch_pos++;
    }
    // Per-patch space, sized to match the maximum size of the patch.
    for (const auto& patch_max_size_pair : data.patch_max_size_) {
      try {
        data.scratch_.push_back(Tensor({""}, {patch_max_size_pair.second}));
      } catch (const std::string& err_msg) {
        throw ERROR_MSG("Failed to call Tensor(). Error:\n\t[", err_msg, "]");
      }
      data.scratch_map_[patch_max_size_pair.first] = patch_pos++;
    }
    // Extra space for storing copies made during cut operations.
    // TODO(martinop): minor optimizations possible: When consecutive cuts apply
    // to the same grid tensor, only one copy needs to be stored.
    for (const auto& cut_copy_size_pair : cut_copy_size) {
      try {
        data.scratch_.push_back(Tensor({""}, {cut_copy_size_pair.second}));
      } catch (const std::string& err_msg) {
        throw ERROR_MSG("Failed to call Tensor(). Error:\n\t[", err_msg, "]");
      }
      data.scratch_map_[cut_copy_size_pair.first] = patch_pos++;
    }
  }

  return data;
}

void ContractionData::ContractGrid(
    std::list<ContractionOperation> ordering, std::size_t output_index,
    std::unordered_map<std::string, bool> active_patches) {
  Tensor* output{nullptr};

  while (!ordering.empty()) {
    const ContractionOperation op = ordering.front();
    ordering.pop_front();
    switch (op.op_type) {
      case ContractionOperation::EXPAND: {
        // Multiply by the new tensor and store result in scratch.
        Tensor& next =
            (*tensor_grid_)[op.expand.tensor[0]][op.expand.tensor[1]];
        if (!active_patches[op.expand.id]) {
          // First tensor in patch.
          get_scratch(op.expand.id) = next;
          active_patches[op.expand.id] = true;
          continue;
        }
        Tensor& prev = get_scratch(op.expand.id);
        std::string result_id = result_space(patch_max_size_[op.expand.id]);
        Tensor& result = get_scratch(result_id);
        try {
          multiply(prev, next, result, get_scratch(kGeneralSpace).data());
        } catch (const std::string& err_msg) {
          throw ERROR_MSG("Failed to call multiply(). Error:\n\t[", err_msg,
                          "]");
        }
        if (ordering.empty()) {
          output = &result;
          continue;
        }
        std::size_t temp = scratch_map_[op.expand.id];
        scratch_map_[op.expand.id] = scratch_map_[result_id];
        scratch_map_[result_id] = temp;
        continue;
      }
      case ContractionOperation::CUT: {
        const std::string index = index_name(op.cut.tensors);
        // If no error is caught, index will be initialized.
        Tensor& tensor_a =
            (*tensor_grid_)[op.cut.tensors[0][0]][op.cut.tensors[0][1]];
        Tensor& copy_a =
            get_scratch(cut_copy_name(op.cut.tensors, /*side = */ 0));
        copy_a = tensor_a;
        // List of values to evaluate on the cut.
        std::vector<std::size_t> values = op.cut.values;
        if (values.empty()) {
          std::size_t num_values = tensor_a.get_index_to_dimension().at(index);
          for (std::size_t i = 0; i < num_values; ++i) {
            values.push_back(i);
          }
        }
        if (op.cut.tensors.size() > 1) {
          // This is a normal cut; each value adds to the same amplitude.
          Tensor& tensor_b =
              (*tensor_grid_)[op.cut.tensors[1][0]][op.cut.tensors[1][1]];
          Tensor& copy_b =
              get_scratch(cut_copy_name(op.cut.tensors, /*side = */ 1));
          copy_b = tensor_b;
          for (std::size_t val : values) {
            try {
              copy_a.project(index, val, tensor_a);
            } catch (const std::string& err_msg) {
              throw ERROR_MSG("Failed to call project(). Error:\n\t[", err_msg,
                              "]");
            }
            try {
              copy_b.project(index, val, tensor_b);
            } catch (const std::string& err_msg) {
              throw ERROR_MSG("Failed to call project(). Error:\n\t[", err_msg,
                              "]");
            }
            ContractGrid(ordering, output_index, active_patches);
          }
          tensor_a = copy_a;
          tensor_b = copy_b;
        } else {
          // This is a terminal cut; each value adds to a different amplitude.
          output_index *= values.size();
          for (std::size_t val : values) {
            try {
              copy_a.project(index, val, tensor_a);
            } catch (const std::string& err_msg) {
              throw ERROR_MSG("Failed to call project(). Error:\n\t[", err_msg,
                              "]");
            }
            ContractGrid(ordering, output_index, active_patches);
            output_index++;
          }
          tensor_a = copy_a;
        }
        return;  // Post-cut contraction is handled recursively.
      }
      case ContractionOperation::MERGE: {
        // Multiply two existing tensors and store result in scratch.
        Tensor& patch_1 = get_scratch(op.merge.source_id);
        Tensor& patch_2 = get_scratch(op.merge.target_id);
        if (!active_patches[op.merge.target_id]) {
          // Copy the old patch into the new space.
          patch_2 = patch_1;
          active_patches[op.merge.target_id] = true;
          continue;
        }
        std::string result_id =
            result_space(patch_max_size_[op.merge.target_id]);
        Tensor& result = get_scratch(result_id);
        try {
          multiply(patch_1, patch_2, result, get_scratch(kGeneralSpace).data());
        } catch (const std::string& err_msg) {
          throw ERROR_MSG("Failed to call multiply(). Error:\n\t[", err_msg,
                          "]");
        }
        if (ordering.empty()) {
          output = &result;
          continue;
        }
        std::size_t temp = scratch_map_[op.merge.target_id];
        scratch_map_[op.merge.target_id] = scratch_map_[result_id];
        scratch_map_[result_id] = temp;
        continue;
      }
    }
  }
  // Check for an output of size larger than 1
  if (output->size() != 1) {
    throw ERROR_MSG("Contraction did not complete; final tensor is: ",
                    output->tensor_to_string());
  }
  (*amplitudes_)[output_index] += (*output->data());
  return;
}

// External methods

void ordering_data_to_contraction_ordering(
    const QflexInput& input, std::list<ContractionOperation>* ordering) {
  if (ordering == nullptr) {
    throw ERROR_MSG("Ordering must be non-null.");
  }

  static const std::regex cut_value_regex("\\([0-9,]*\\)");
  std::string line;
  std::string operation;

  auto error_msg = [&line](const auto& msg) {
    return ERROR_MSG("Parsing failed on line: '", line, "' with error: ", msg);
  };

  // Lambda function to check if indexes are in the right range
  auto check_index = [&input, &error_msg](auto index) {
    using index_type =
        std::remove_reference_t<std::remove_cv_t<decltype(index)>>;
    const std::size_t grid_size = input.grid.I * input.grid.J;
    if (static_cast<std::size_t>(static_cast<index_type>(grid_size)) !=
        grid_size) {
      throw error_msg("Grid is too big.");
    } else if (index < 0 or index >= static_cast<index_type>(grid_size)) {
      throw error_msg(
          utils::concat("Index ", index, " must be within grid boundaries."));
    }
  };

  // Read one line at a time from the ordering, skipping comments.
  for (const auto& line : input.ordering.instructions) {
    if (line.empty() || line[0] == '#') continue;
    std::stringstream ss(line);
    // The first element is the operation (expand, cut, or merge).
    ss >> operation;
    if (operation == "expand") {
      std::string patch;
      long int index;
      ss >> patch;
      ss >> index;

      // Check index
      check_index(index);

      std::vector<std::size_t> position = get_qubit_coords(index, input.grid.J);
      if (find_grid_coord_in_list(input.grid.qubits_off, position[0],
                                  position[1])) {
        throw error_msg("Index must specify an active qubit.");
      }
      ordering->emplace_back(ExpandPatch(patch, position));

    } else if (operation == "cut") {
      std::vector<std::size_t> values;
      std::string values_str;
      std::smatch match;
      ss >> values_str;
      if (!std::regex_match(values_str, match, cut_value_regex)) {
        throw error_msg(
            "Cut values must be comma-separated ints, e.g. (0,1,3).");
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
      long int index_1, index_2;
      ss >> index_1;

      // Check index
      check_index(index_1);

      std::vector<std::size_t> position_1 =
          get_qubit_coords(index_1, input.grid.J);
      if (find_grid_coord_in_list(input.grid.qubits_off, position_1[0],
                                  position_1[1])) {
        throw error_msg("Index 1 must specify an active qubit.");
      }
      if (ss.eof()) {
        ordering->emplace_back(CutIndex({position_1}, values));
      } else {
        ss >> index_2;

        // Check index
        check_index(index_2);

        std::vector<std::size_t> position_2 =
            get_qubit_coords(index_2, input.grid.J);
        if (find_grid_coord_in_list(input.grid.qubits_off, position_2[0],
                                    position_2[1])) {
          throw error_msg("Index 2 must specify an active qubit.");
        }
        // If indices are listed in reverse order, swap them to prevent issues
        // in tensor contraction.
        if (index_1 < index_2) {
          ordering->emplace_back(CutIndex({position_1, position_2}, values));
        } else {
          ordering->emplace_back(CutIndex({position_2, position_1}, values));
        }
      }
    } else if (operation == "merge") {
      std::string patch_1, patch_2;
      ss >> patch_1;
      ss >> patch_2;
      ordering->emplace_back(MergePatches(patch_1, patch_2));
    } else {
      throw error_msg("Received an invalid operation in config.");
    }
  }

  // Ensure ordering generated is valid
  ValidateOrdering(*ordering);
}

std::string index_name(const std::vector<std::size_t>& p1,
                       const std::vector<std::size_t>& p2) {
  if (p1.size() == 2 && p2.size() == 2) {
    // Two-qubit contraction.
    return utils::concat("(", p1[0], ",", p1[1], "),(", p2[0], ",", p2[1], ")");
  }
  if (p1.size() == 3 && p2.size() == 3) {
    // Single-qubit contraction, or virtual index.
    return utils::concat("(", p1[0], ",", p1[1], ",", p1[2], "),(", p2[0], ",",
                         p2[1], ",", p2[2], ")");
  }
  if (p1.size() == 2 && p2.empty()) {
    // Final qubit output value assignment.
    return utils::concat("(", p1[0], ",", p1[1], "),(o)");
  }

  std::stringstream ss;
  ss << "Failed to construct tensor name with the following vectors: p1 = [";
  for (const auto& p : p1) ss << std::to_string(p) << ',';
  ss << "] and p2 = [";
  for (const auto& p : p2) ss << std::to_string(p) << ',';
  ss << "].";
  throw ERROR_MSG(ss.str());
}

std::string index_name(const std::vector<std::vector<std::size_t>>& tensors) {
  if (tensors.size() == 2) {
    try {
      return index_name(tensors.at(0), tensors.at(1));
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call index_name(). Error:\n\t[", err_msg, "]");
    }
  }
  if (tensors.size() == 1) {
    try {
      return index_name(tensors.at(0), {});
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call index_name(). Error:\n\t[", err_msg, "]");
    }
  }
  throw ERROR_MSG("Failed to construct tensor name with input tensors size: ",
                  tensors.size(), ".");
}

std::vector<std::size_t> get_qubit_coords(std::size_t q, std::size_t J) {
  std::size_t i = q / J;
  return std::vector<std::size_t>({i, q - i * J});
}

bool find_grid_coord_in_list(
    const std::optional<std::vector<std::vector<std::size_t>>>& coord_list,
    const std::size_t i, const std::size_t j) {
  return coord_list.has_value() &&
         find(coord_list.value().begin(), coord_list.value().end(),
              std::vector<std::size_t>({i, j})) != coord_list.value().end();
}

void ValidateOrdering(const std::list<ContractionOperation>& ordering) {
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
  std::string error_msg;
  for (const auto& op : ordering) {
    switch (op.op_type) {
      case ContractionOperation::EXPAND: {
        if (patches[op.expand.id].is_used)
          error_msg = utils::concat(
              error_msg, "\nTensor at (", op.expand.tensor[0], ",",
              op.expand.tensor[1], ") is added to non-empty patch ",
              op.expand.id.c_str(), " after a cut.");
        if (patches[op.expand.id].is_merged)
          error_msg = utils::concat(
              error_msg, "\nTensor at (", op.expand.tensor[0], ",",
              op.expand.tensor[1], ") is added to non-empty patch ",
              op.expand.id.c_str(), " after a cut.");

        std::string tensor_name = utils::concat("(", op.expand.tensor[0], ",",
                                                op.expand.tensor[1], ")");
        if (used_tensors.find(tensor_name) != used_tensors.end())
          error_msg = utils::concat(error_msg, "\nTensor ", tensor_name,
                                    " is contracted multiple times.");

        used_tensors.insert(tensor_name);
        patches[op.expand.id].is_active = true;
        continue;
      }
      case ContractionOperation::CUT: {
        try {
          for (auto& patch_pair : patches) {
            if (patch_pair.second.is_active) {
              patch_pair.second.is_used = true;
            }
          }
          const std::string index = index_name(op.cut.tensors);
          // If no error is caught, index will be initialized.
          if (cut_indices.find(index) != cut_indices.end())
            error_msg = utils::concat(error_msg, "\nIndex ", index.c_str(),
                                      " is cut multiple times.");

          cut_indices.insert(index);
        } catch (const std::string& err_msg) {
          throw ERROR_MSG("Failed to during CUT. Error:\n\t[", err_msg, "]");
        }
        continue;
      }
      case ContractionOperation::MERGE: {
        if (patches[op.merge.source_id].is_merged)
          error_msg =
              utils::concat(error_msg, "\nPatch ", op.merge.source_id.c_str(),
                            " is merged multiple times.");

        if (patches[op.merge.target_id].is_used)
          error_msg =
              utils::concat(error_msg, "\nPatch ", op.merge.source_id.c_str(),
                            " is merged into non-empty patch ",
                            op.merge.target_id.c_str(), " after a cut.");

        patches[op.merge.source_id].is_merged = true;
        patches[op.merge.target_id].is_active = true;
        continue;
      }
    }
  }

  if (not std::empty(error_msg)) throw error_msg;
}

void ContractGrid(const std::list<ContractionOperation>& ordering,
                  std::vector<std::vector<Tensor>>* tensor_grid,
                  std::vector<std::complex<double>>* amplitudes) {
  if (tensor_grid == nullptr) {
    throw ERROR_MSG("Tensor grid must be non-null.");
  }
  if (amplitudes == nullptr) {
    throw ERROR_MSG("Amplitude return vector must be non-null.");
  }

  // Populate ContractionData and perform grid contraction.
  ContractionData data;
  try {
    data = ContractionData::Initialize(ordering, tensor_grid, amplitudes);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call Initialize(). Error:\n\t[", err_msg, "]");
  }
  // If no error is caught, data will be initialized.
  std::unordered_map<std::string, bool> active_patches;
  for (const auto& patch : data.scratch_list()) {
    active_patches[patch] = false;
  }
  try {
    data.ContractGrid(ordering, /*output_index = */ 0, active_patches);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call ContractGrid(). Error:\n\t[", err_msg, "]");
  }
}

}  // namespace qflex
