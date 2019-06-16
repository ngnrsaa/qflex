#include "contraction_utils.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <unordered_set>

constexpr char kResultSpace[] = "null";

namespace {

// Copies a ContractionOrdering for reuse.
ContractionOrdering copy_order(const ContractionOrdering& ordering) {
  ContractionOrdering new_order;
  for (const auto& op : ordering) {
    switch (op->op_type) {
      case ContractionOperation::EXPAND: {
        const auto* expand = dynamic_cast<const ExpandPatch*>(op.get());
        new_order.emplace_back(new ExpandPatch(*expand));
        break;
      }
      case ContractionOperation::CUT: {
        const auto* cut = dynamic_cast<const CutIndex*>(op.get());
        new_order.emplace_back(new CutIndex(*cut));
        break;
      }
      case ContractionOperation::MERGE: {
        const auto* merge = dynamic_cast<const MergePatches*>(op.get());
        new_order.emplace_back(new MergePatches(*merge));
        break;
      }
    }
  }
  return new_order;
}

}  // namespace

void ContractionData::ContractGrid(
    ContractionOrdering ordering, int output_index,
    std::unordered_map<std::string, bool> active_patches) {
  while (!ordering.empty()) {
    std::unique_ptr<ContractionOperation> op = std::move(ordering.front());
    ordering.pop_front();
    switch (op->op_type) {
      case ContractionOperation::EXPAND: {
        const auto* expand = dynamic_cast<const ExpandPatch*>(op.get());
        // Multiply by the new tensor and store result in scratch.
        MKLTensor& next = (*tensor_grid_)[expand->tensor[0]][expand->tensor[1]];
        if (!active_patches[expand->id]) {
          // First tensor in patch.
          scratch_[scratch_map_[expand->id]] = next;
          active_patches[expand->id] = true;
          continue;
        }
        MKLTensor& prev = scratch_[scratch_map_[expand->id]];
        MKLTensor& result = scratch_[scratch_map_[kResultSpace]];
        multiply(prev, next, result, scratch_[0].data());
        if (ordering.empty()) continue;
        int temp = scratch_map_[expand->id];
        scratch_map_[expand->id] = scratch_map_[kResultSpace];
        scratch_map_[kResultSpace] = temp;
        continue;
      }
      case ContractionOperation::CUT: {
        // DO NOT SUBMIT: allocating copies may be overly-costly!
        // (add extra scratch?)
        const auto* cut = dynamic_cast<const CutIndex*>(op.get());
        const std::string index = index_name(cut->tensors);
        MKLTensor* tensor_a =
            &(*tensor_grid_)[cut->tensors[0][0]][cut->tensors[0][1]];
        const MKLTensor copy_a(*tensor_a);
        // List of values to evaluate on the cut.
        std::vector<int> values = cut->values;
        if (values.empty()) {
          int num_values = tensor_a->get_index_to_dimension().at(index);
          for (int i = 0; i < num_values; ++ i) {
            values.push_back(i);
          }
        }
        if (cut->tensors.size() > 1) {
          // This is a normal cut; each value adds to the same amplitude.
          MKLTensor* tensor_b =
              &(*tensor_grid_)[cut->tensors[1][0]][cut->tensors[1][1]];
          const MKLTensor copy_b(*tensor_b);
          for (int val : values) {
            copy_a.project(index_name(cut->tensors), val, *tensor_a);
            copy_b.project(index_name(cut->tensors), val, *tensor_b);
            ContractGrid(copy_order(ordering), output_index, active_patches);
          }
          *tensor_a = copy_a;
          *tensor_b = copy_b;
        } else {
          // This is a terminal cut; each value adds to a different amplitude.
          output_index *= values.size();
          for (int val : values) {
            copy_a.project(index_name(cut->tensors), val, *tensor_a);
            ContractGrid(copy_order(ordering), output_index, active_patches);
            output_index++;
          }
          *tensor_a = copy_a;
        }
        return;  // Post-cut contraction is handled recursively.
      }
      case ContractionOperation::MERGE: {
        const auto* merge = dynamic_cast<const MergePatches*>(op.get());
        // Multiply two existing tensors and store result in scratch.
        MKLTensor& patch_1 = scratch_[scratch_map_[merge->source_id]];
        MKLTensor& patch_2 = scratch_[scratch_map_[merge->target_id]];
        if (!active_patches[merge->target_id]) {
          // Copy the old patch into the new space.
          patch_2 = patch_1;
          active_patches[merge->target_id] = true;
          continue;
        }
        MKLTensor& result = scratch_[scratch_map_[kResultSpace]];
        multiply(patch_1, patch_2, result, scratch_[0].data());
        if (ordering.empty()) continue;
        int temp = scratch_map_[merge->target_id];
        scratch_map_[merge->target_id] = scratch_map_[kResultSpace];
        scratch_map_[kResultSpace] = temp;
        continue;
      }
    }
  }
  MKLTensor& result = scratch_[scratch_map_[kResultSpace]];
  if (result.size() != 1) {
    std::cout << "Contraction did not complete; final tensor is ";
    result.print();
    assert(false);
  }
  (*amplitudes_)[output_index] += (*result.data());
  return;
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
  assert(false && "Failed to construct tensor name.");
  return "";
}

std::string index_name(const std::vector<std::vector<int>>& tensors) {
  if (tensors.size() == 2) {
    return index_name(tensors.at(0), tensors.at(1));
  }
  if (tensors.size() == 1) {
    return index_name(tensors.at(0), {});
  }
  assert(false && "Failed to construct tensor name.");
  return "";
}

bool IsOrderingValid(const ContractionOrdering& ordering) {
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
    switch (op->op_type) {
      case ContractionOperation::EXPAND: {
        const auto* expand = dynamic_cast<const ExpandPatch*>(op.get());
        if (patches[expand->id].is_used) {
          next_error += snprintf(
              error_msg + next_error, error_space - next_error,
              "Tensor at (%d,%d) is added to non-empty patch %s after a cut.\n",
              expand->tensor[0], expand->tensor[1], expand->id.c_str());
        }
        if (patches[expand->id].is_merged) {
          next_error += snprintf(
              error_msg + next_error, error_space - next_error,
              "Tensor at (%d,%d) is added to previously-merged patch %s.\n",
              expand->tensor[0], expand->tensor[1], expand->id.c_str());
        }
        char tensor_name[20];
        snprintf(tensor_name, sizeof(tensor_name), "(%d,%d)", expand->tensor[0],
                 expand->tensor[1]);
        if (used_tensors.find(tensor_name) != used_tensors.end()) {
          next_error += snprintf(
              error_msg + next_error, error_space - next_error,
              "Tensor %s is contracted multiple times.\n", tensor_name);
        }
        used_tensors.insert(tensor_name);
        patches[expand->id].is_active = true;
        continue;
      }
      case ContractionOperation::CUT: {
        const auto* cut = dynamic_cast<const CutIndex*>(op.get());
        for (auto& patch_pair : patches) {
          if (patch_pair.second.is_active) {
            patch_pair.second.is_used = true;
          }
        }
        const std::string index = index_name(cut->tensors);
        if (cut_indices.find(index) != cut_indices.end()) {
          next_error +=
              snprintf(error_msg + next_error, error_space - next_error,
                       "Index %s is cut multiple times.\n", index.c_str());
        }
        cut_indices.insert(index);
        continue;
      }
      case ContractionOperation::MERGE: {
        const auto* merge = dynamic_cast<const MergePatches*>(op.get());
        if (patches[merge->source_id].is_merged) {
          next_error += snprintf(
              error_msg + next_error, error_space - next_error,
              "Patch %s is merged multiple times.\n", merge->source_id.c_str());
        }
        if (patches[merge->target_id].is_used) {
          next_error += snprintf(
              error_msg + next_error, error_space - next_error,
              "Patch %s is merged into non-empty patch %s after a cut.\n",
              merge->source_id.c_str(), merge->target_id.c_str());
        }
        patches[merge->source_id].is_merged = true;
        patches[merge->target_id].is_active = true;
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

void ContractGrid(const ContractionOrdering& ordering, const int data_size,
                  std::vector<std::vector<MKLTensor>>* tensor_grid,
                  std::vector<std::complex<double>>* amplitudes) {
  assert(amplitudes != nullptr && "Amplitude return vector must be non-null.");
  assert(IsOrderingValid(ordering));

  // Acquire the necessary scratch space and populate ContractionData.
  std::unordered_map<std::string, bool> active_patches;
  ContractionData data;
  data.max_rank_ = data_size;
  data.tensor_grid_ = tensor_grid;
  data.amplitudes_ = amplitudes;
  // Scratch space for reordering during multiplication.
  data.scratch_.push_back(MKLTensor({""}, {data.max_rank_}));

  // Scratch space for temporarily storing multiplication results.
  data.scratch_.push_back(MKLTensor({""}, {data.max_rank_}));
  data.scratch_map_[kResultSpace] = 1;

  // Scratch space for storing patch tensors.
  int i = data.scratch_.size();
  for (const auto& op : ordering) {
    std::string mutate_id;
    if (op->op_type == ContractionOperation::EXPAND) {
      const auto* expand = dynamic_cast<const ExpandPatch*>(op.get());
      mutate_id = expand->id;
    } else if (op->op_type == ContractionOperation::MERGE) {
      const auto* merge = dynamic_cast<const MergePatches*>(op.get());
      mutate_id = merge->target_id;
    } else {
      continue;
    }
    if (data.scratch_map_.find(mutate_id) == data.scratch_map_.end()) {
      data.scratch_map_[mutate_id] = i++;
      active_patches[mutate_id] = false;
      // TODO(martinop): reduce scratch space usage when possible.
      data.scratch_.push_back(MKLTensor({""}, {data.max_rank_}));
    }
  }
  data.ContractGrid(copy_order(ordering), /*output_index = */ 0,
                    active_patches);
}
