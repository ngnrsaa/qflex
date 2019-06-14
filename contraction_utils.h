#include <list>
#ifndef CONTRACTION_UTILS_
#define CONTRACTION_UTILS_

#include <complex>
#include <memory>
#include <vector>

#include "mkl_tensor.h"

struct ContractionOperation {
 public:
  enum OpType { EXPAND, CUT, MERGE };
  OpType op_type;
  virtual ~ContractionOperation() {}

 protected:
  ContractionOperation(OpType op_type) : op_type(op_type) {}
};

struct ExpandPatch : public ContractionOperation {
  ExpandPatch(char id, std::vector<int> tensor)
      : ContractionOperation(EXPAND), id(id), tensor(tensor) {}
  ~ExpandPatch() override {}
  // ID of the patch expanded by this operation.
  char id;
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
  MergePatches(char source_id, char target_id)
      : ContractionOperation(MERGE),
        source_id(source_id),
        target_id(target_id) {}
  ~MergePatches() override {}
  // IDs of the patches being merged. The resulting patch will use target_id.
  char source_id;
  char target_id;
};

using ContractionOrdering = std::list<std::unique_ptr<ContractionOperation>>;

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
 * */
bool IsOrderingValid(const ContractionOrdering& ordering);

#endif  // CONTRACTION_UTILS_
