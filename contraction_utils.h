#ifndef CONTRACTION_UTILS_
#define CONTRACTION_UTILS_

#include <memory>
#include <vector>

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
  char id;
  std::vector<int> tensor;
};

struct CutIndex : public ContractionOperation {
  CutIndex(std::vector<std::vector<int>> tensors)
      : ContractionOperation(CUT), tensors(tensors) {}
  CutIndex(std::vector<int> tensor_a, std::vector<int> tensor_b)
      : ContractionOperation(CUT), tensors({tensor_a, tensor_b}) {}
  ~CutIndex() override {}
  std::vector<std::vector<int>> tensors;
};

struct MergePatches : public ContractionOperation {
  MergePatches(char source_id, char target_id)
      : ContractionOperation(MERGE),
        source_id(source_id),
        target_id(target_id) {}
  ~MergePatches() override {}
  char source_id;
  char target_id;
};

using ContractionOrdering = std::vector<std::unique_ptr<ContractionOperation>>;

#endif  // CONTRACTION_UTILS_
