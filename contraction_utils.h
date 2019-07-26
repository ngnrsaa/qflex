#ifndef CONTRACTION_UTILS_
#define CONTRACTION_UTILS_

#include <complex>
#include <fstream>
#include <iostream>
#include <list>
#include <memory>
#include <vector>

#include "mkl_tensor.h"

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
std::string index_name(const std::vector<int>& p1, const std::vector<int>& p2);

// As above, but accepts a list of qubit locations. If the list has only one
// element, an empty vector will be provided as the second element.
std::string index_name(const std::vector<std::vector<int>>& tensors);

/**
 * Returns spatial coordinates on the grid for the given qubit q.
 * @param q int with the qubit number.
 * @param J int with the second spatial dimension of the grid of qubits.
 * @return vector<int> with the spatial coordinates of qubit q on the grid.
 */
std::vector<int> get_qubit_coords(int q, int J);

/**
 * Helper function to find a grid coordinate in a list of coordinates.
 * @param coord_list optional vector of qubit positions.
 * @param i int with the first spatial dimension of the target position.
 * @param j int with the second spatial dimension of the target position.
 * @return true if (i, j) is in coord_list.
 */
bool find_grid_coord_in_list(
    const std::optional<std::vector<std::vector<int>>>& coord_list, const int i,
    const int j);

struct ContractionOperation {
 public:
  enum OpType { EXPAND, CUT, MERGE };
  OpType op_type;
  virtual ~ContractionOperation() {}

 protected:
  ContractionOperation(OpType op_type) : op_type(op_type) {}
};

struct ExpandPatch : public ContractionOperation {
  ExpandPatch(std::string id, std::vector<int> tensor)
      : ContractionOperation(EXPAND), id(id), tensor(tensor) {}
  ~ExpandPatch() override {}
  // ID of the patch expanded by this operation.
  std::string id;
  // Tensor to contract into the patch.
  std::vector<int> tensor;
};

struct CutIndex : public ContractionOperation {
  CutIndex(std::vector<std::vector<int>> tensors, std::vector<int> values = {})
      : ContractionOperation(CUT), tensors(tensors), values(values) {}
  ~CutIndex() override {}
  // Tensors connected by the cut index. If only one tensor is provided, this
  // represents a 'terminal' cut - i.e., a choice of output value for a qubit.
  std::vector<std::vector<int>> tensors;
  // Values to assign to this index when cutting. If left blank, all possible
  // values will be assigned.
  std::vector<int> values;
};

struct MergePatches : public ContractionOperation {
  MergePatches(std::string source_id, std::string target_id)
      : ContractionOperation(MERGE),
        source_id(source_id),
        target_id(target_id) {}
  ~MergePatches() override {}
  // IDs of the patches being merged. The resulting patch will use target_id.
  std::string source_id;
  std::string target_id;
};

using ContractionOrdering = std::list<std::unique_ptr<ContractionOperation>>;

// Copies a ContractionOrdering for reuse.
ContractionOrdering copy_order(const ContractionOrdering& ordering);

/**
 * Parses a grid contraction ordering from the given stream.
 * @param circuit_data std::istream containing ordering as a string.
 * @param I int with the first spatial dimension of the grid of qubits.
 * @param J int with the second spatial dimension of the grid of qubits.
 * @param off vector<vector<int>> with the coords. of the qubits turned off.
 * @param ordering pointer to ContractionOrdering output object.
 * @return false if parsing failed at any point, true otherwise.
 **/
bool ordering_data_to_contraction_ordering(
    std::istream* ordering_data, const int I, const int J,
    const std::optional<std::vector<std::vector<int>>>& off,
    ContractionOrdering* ordering);

// Helper class for the external ContractGrid method. This should not be
// initialized by external users.
class ContractionData {
 public:
  /**
   * Generates a ContractionData object from a contraction ordering, and logs
   * total memory allocated for the contraction.
   * @param ordering ContractionOrdering listing operations to perform.
   * @param tensor_grid 2D vector<MKLTensor> holding the tensor grid. Consumes
   * output from grid_of_tensors_3D_to_2D.
   * @param amplitudes vector of amplitudes for each final output requested.
   */
  static ContractionData Initialize(
      const ContractionOrdering& ordering,
      std::vector<std::vector<MKLTensor>>* tensor_grid,
      std::vector<std::complex<double>>* amplitudes);

  // Keys for scratch-space tensors. Do not reuse outside this file.
  static constexpr char kGeneralSpace[] = "_general_internal_";
  static constexpr char kResultSpace[] = "_result_internal_";

  /**
   * Recursive helper for the external ContractGrid method below. This method
   * calls itself recursively on each "cut" operation.
   * @param ordering ContractionOrdering listing operations to perform.
   * @param output_index int marking which amplitude will be updated next.
   * @param active_patches list of patches already created in scratch space.
   */
  void ContractGrid(ContractionOrdering ordering, int output_index,
                    std::unordered_map<std::string, bool> active_patches);

  MKLTensor& get_scratch(std::string id) { return scratch_[scratch_map_[id]]; }

  /**
   * Assigns a name to a given cut-copy tensor for mapping into scratch space.
   * @param index the index being cut.
   * @param side the side of the cut referenced.
   * @return the key used for this cut-copy.
   */
  static std::string cut_copy_name(std::vector<std::vector<int>> index,
                                   int side) {
    std::string base = index_name(index);
    char buffer[64];
    int len =
        snprintf(buffer, sizeof(buffer), "cut-%s:side-%d", base.c_str(), side);
    return std::string(buffer, len);
  }

  // Gets the index of result scratch space of the given rank.
  static std::string result_space(int rank) {
    char buffer[64];
    int len = snprintf(buffer, sizeof(buffer), "%s%d", kResultSpace, rank);
    return std::string(buffer, len);
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
  std::vector<MKLTensor> scratch_;

  // Map of patch IDs/index names to scratch locations.
  std::unordered_map<std::string, int> scratch_map_;

  // Max rank of each patch generated during contraction.
  std::unordered_map<std::string, int> patch_rank_;

  // Contains the tensor grid produced by grid_of_tensors_3D_to_2D.
  std::vector<std::vector<MKLTensor>>* tensor_grid_;

  // Amplitudes for each final output requested.
  std::vector<std::complex<double>>* amplitudes_;
};

/**
 * Performs basic sanity checks on the given contraction ordering:
 *   - A patch cannot be expanded after being merged into another patch.
 *     - Workaround: expand the target patch instead.
 *   - A patch cannot be expanded/merged into both before and after a cut.
 *     This includes ALL cuts, not just ones adjacent to the patch.
 *     - Workaround: use a merge to rename the patch after each cut.
 *   - Each tensor can only be in one ExpandPatch.
 *   - Each index can only be in one CutIndex.
 *   - Each patch can only be the source in one MergePatches.
 * @param ordering ContractionOrdering listing operations to perform.
 * @return true is the ordering is valid, false otherwise.
 */
bool IsOrderingValid(const ContractionOrdering& ordering);

/**
 * Performs contraction operations specified by 'ordering' on tensor_grid.
 *
 * This method will allocate (data_size)*(# of patches + 2) units of memory
 * for performing the grid contraction. A small amout of additional space is
 * allocated during cuts to preserve tensor_grid values.
 * @param ordering ContractionOrdering listing operations to perform.
 * @param tensor_grid 2D vector<MKLTensor> holding the tensor grid. Consumes
 * output from grid_of_tensors_3D_to_2D.
 * @param amplitudes vector of amplitudes for each final output requested.
 */
void ContractGrid(const ContractionOrdering& ordering,
                  std::vector<std::vector<MKLTensor>>* tensor_grid,
                  std::vector<std::complex<double>>* amplitudes);

}  // namespace qflex

#endif  // CONTRACTION_UTILS_
