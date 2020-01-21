#ifndef CONTRACTION_UTILS_
#define CONTRACTION_UTILS_

#include <complex>
#include <fstream>
#include <iostream>
#include <list>
#include <memory>
#include <vector>

#include "global.h"
#include "input.h"
#include "tensor.h"

namespace qflex {

/**
 * Method for assigning names to indices in the grid. Accepted formats include:
 *   - {i_1, j_1}, {i_2, j_2}: two-qubit contraction.
 *   - {i_1, j_1, k_1}, {i_2, j_2, k_2}: single-qubit gate or virtual index.
 *   - {i, j}, {}: output-value assignment.
 * @param p1 position of the first connected tensor.
 * @param p2 position of the second connected tensor.
 * @return string name of the index.
 */
std::string index_name(const std::vector<std::size_t>& p1,
                       const std::vector<std::size_t>& p2);

// As above, but accepts a list of qubit locations. If the list has only one
// element, an empty vector will be provided as the second element.
std::string index_name(const std::vector<std::vector<std::size_t>>& tensors);

/**
 * Returns spatial coordinates on the grid for the given qubit q.
 * @param q std::size_t with the qubit number.
 * @param J std::size_t with the second spatial dimension of the grid of qubits.
 * @return vector<std::size_t> with the spatial coordinates of qubit q on the
 * grid.
 */
std::vector<std::size_t> get_qubit_coords(std::size_t q, std::size_t J);

/**
 * Helper function to find a grid coordinate in a list of coordinates.
 * @param coord_list optional vector of qubit positions.
 * @param i std::size_t with the first spatial dimension of the target position.
 * @param j std::size_t with the second spatial dimension of the target
 * position.
 * @return true if (i, j) is in coord_list.
 */
bool find_grid_coord_in_list(
    const std::optional<std::vector<std::vector<std::size_t>>>& coord_list,
    const std::size_t i, const std::size_t j);

struct ExpandPatch {
  ExpandPatch() {}
  ExpandPatch(std::string id, std::vector<std::size_t> tensor)
      : id(id), tensor(tensor) {}
  // ID of the patch expanded by this operation.
  std::string id;
  // Tensor to contract into the patch.
  std::vector<std::size_t> tensor;
};

struct CutIndex {
  CutIndex() {}
  CutIndex(std::vector<std::vector<std::size_t>> tensors,
           std::vector<std::size_t> values = {})
      : tensors(tensors), values(values) {}
  // Tensors connected by the cut index. If only one tensor is provided, this
  // represents a 'terminal' cut - i.e., a choice of output value for a qubit.
  std::vector<std::vector<std::size_t>> tensors;
  // Values to assign to this index when cutting. If left blank, all possible
  // values will be assigned.
  std::vector<std::size_t> values;
};

struct MergePatches {
  MergePatches() {}
  MergePatches(std::string source_id, std::string target_id)
      : source_id(source_id), target_id(target_id) {}
  // IDs of the patches being merged. The resulting patch will use target_id.
  std::string source_id;
  std::string target_id;
};

struct ContractionOperation {
 public:
  ContractionOperation(ExpandPatch expand) : op_type(EXPAND), expand(expand) {}
  ContractionOperation(CutIndex cut) : op_type(CUT), cut(cut) {}
  ContractionOperation(MergePatches merge) : op_type(MERGE), merge(merge) {}

  enum OpType { EXPAND, CUT, MERGE };
  OpType op_type;

  // Only one of these should be populated, matching op_type.
  const ExpandPatch expand;
  const CutIndex cut;
  const MergePatches merge;
};

/**
 * Parses a grid contraction ordering from the given stream.
 * @param QflexInput containing simulation information
 * @param ordering pointer to std::list<ContractionOperation> output object.
 **/
void ordering_data_to_contraction_ordering(
    const QflexInput& input, std::list<ContractionOperation>* ordering);

// Helper class for the external ContractGrid method. This should not be
// initialized by external users.
class ContractionData {
 public:
  /**
   * Generates a ContractionData object from a contraction ordering, and logs
   * total memory allocated for the contraction.
   * @param ordering std::list<ContractionOperation> listing operations to
   * perform.
   * @param tensor_grid 2D vector<Tensor> holding the tensor grid. Consumes
   * output from grid_of_tensors_3D_to_2D.
   * @param amplitudes vector of amplitudes for each final output requested.
   */
  static ContractionData Initialize(
      const std::list<ContractionOperation>& ordering,
      std::vector<std::vector<Tensor>>* tensor_grid,
      std::vector<std::complex<double>>* amplitudes);

  // Keys for scratch-space tensors. Do not reuse outside this file.
  static constexpr char kGeneralSpace[] = "_general_internal_";
  static constexpr char kResultSpace[] = "_result_internal_";

  /**
   * Recursive helper for the external ContractGrid method below. This method
   * calls itself recursively on each "cut" operation.
   * @param ordering std::list<ContractionOperation> listing operations to
   * perform.
   * @param output_index std::size_t marking which amplitude will be updated
   * next.
   * @param active_patches list of patches already created in scratch space.
   */
  void ContractGrid(std::list<ContractionOperation> ordering,
                    std::size_t output_index,
                    std::unordered_map<std::string, bool> active_patches);

  // Gets a reference to the tensor in scratch with ID id.
  Tensor& get_scratch(std::string id) { return scratch_[scratch_map_[id]]; }

  /**
   * Assigns a name to a given cut-copy tensor for mapping into scratch space.
   * @param index the index being cut.
   * @param side the side of the cut referenced.
   * @return the key used for this cut-copy.
   */
  static std::string cut_copy_name(std::vector<std::vector<std::size_t>> index,
                                   std::size_t side) {
    return utils::concat("cut-", index_name(index), ":size-", side);
  }

  // Gets the index of result scratch space of the given size.
  static std::string result_space(std::size_t size) {
    return utils::concat(kResultSpace, size);
  }

  // Gets the names of all scratch tensors.
  std::vector<std::string> scratch_list() {
    std::vector<std::string> names;
    for (const auto& name_pos_pair : scratch_map_) {
      names.push_back(name_pos_pair.first);
    }
    return names;
  }

 private:
  // List of tensors used for scratch space or storage between recursive calls.
  std::vector<Tensor> scratch_;

  // Map of patch IDs/index names to scratch locations on scratch_.
  std::unordered_map<std::string, std::size_t> scratch_map_;

  // Max size of each patch generated during contraction.
  std::unordered_map<std::string, std::size_t> patch_max_size_;

  // Contains the tensor grid produced by grid_of_tensors_3D_to_2D.
  std::vector<std::vector<Tensor>>* tensor_grid_;

  // Amplitudes for each final output requested.
  std::vector<std::complex<double>>* amplitudes_;
};

/**
 * Performs basic sanity checks on the given contraction ordering and throws in
 * case of an error:
 *   - A patch cannot be expanded after being merged into another patch.
 *     - Workaround: expand the target patch instead.
 *   - A patch cannot be expanded/merged into both before and after a cut.
 *     This includes ALL cuts, not just ones adjacent to the patch.
 *     - Workaround: use a merge to rename the patch after each cut.
 *   - Each tensor can only be in one ExpandPatch.
 *   - Each index can only be in one CutIndex.
 *   - Each patch can only be the source in one MergePatches.
 * @param ordering std::list<ContractionOperation> listing operations to
 * perform.
 * @return void.
 */
void ValidateOrdering(const std::list<ContractionOperation>& ordering);

/**
 * Performs contraction operations specified by 'ordering' on tensor_grid.
 *
 * This method will allocate (data_size)*(# of patches + 2) units of memory
 * for performing the grid contraction. A small amout of additional space is
 * allocated during cuts to preserve tensor_grid values.
 * @param ordering std::list<ContractionOperation> listing operations to
 * perform.
 * @param tensor_grid 2D vector<Tensor> holding the tensor grid. Consumes
 * output from grid_of_tensors_3D_to_2D.
 * @param amplitudes vector of amplitudes for each final output requested.
 */
void ContractGrid(const std::list<ContractionOperation>& ordering,
                  std::vector<std::vector<Tensor>>* tensor_grid,
                  std::vector<std::complex<double>>* amplitudes);

}  // namespace qflex

#endif  // CONTRACTION_UTILS_
