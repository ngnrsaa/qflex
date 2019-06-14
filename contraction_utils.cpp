#include "contraction_utils.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <unordered_set>

constexpr char kResultSpace = '\0';

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
  std::unordered_map<char, PatchState> patches;
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
              "Tensor at (%d,%d) is added to non-empty patch %c after a cut.\n",
              expand->tensor[0], expand->tensor[1], expand->id);
        }
        if (patches[expand->id].is_merged) {
          next_error += snprintf(
              error_msg + next_error, error_space - next_error,
              "Tensor at (%d,%d) is added to previously-merged patch %c.\n",
              expand->tensor[0], expand->tensor[1], expand->id);
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
              "Patch %c is merged multiple times.\n", merge->source_id);
        }
        if (patches[merge->target_id].is_used) {
          next_error += snprintf(
              error_msg + next_error, error_space - next_error,
              "Patch %c is merged into non-empty patch %c after a cut.\n",
              merge->source_id, merge->target_id);
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
