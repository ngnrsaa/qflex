#include "contraction_utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace qflex {
namespace {

TEST(ContractionTest, IndexNaming) {
  // Standard two-qubit index.
  std::vector<std::vector<int>> index = {{1, 2}, {3, 4}};
  EXPECT_EQ(index_name(index), "(1,2),(3,4)");

  // Virtual index.
  index = {{1, 2, 3}, {4, 5, 6}};
  EXPECT_EQ(index_name(index), "(1,2,3),(4,5,6)");

  // Terminal index.
  index = {{1, 2}};
  EXPECT_EQ(index_name(index), "(1,2),(o)");
}

TEST(ContractionDeathTest, IndexNamingFailures) {
  // Empty input.
  std::vector<std::vector<int>> index = {};
  EXPECT_DEATH(index_name(index), "");

  // Size mismatch.
  index = {{1, 2}, {3, 4, 5}};
  EXPECT_DEATH(index_name(index), "");

  // Invalid size.
  index = {{1, 2, 3, 4}, {5, 6, 7, 8}};
  EXPECT_DEATH(index_name(index), "");

  // Wrong number of tensor positions.
  index = {{1, 2}, {3, 4}, {5, 6}};
  EXPECT_DEATH(index_name(index), "");

  // Wrong-size terminal index.
  index = {{1, 2, 3}};
  EXPECT_DEATH(index_name(index), "");
}

TEST(ContractionTest, OperationHandling) {
  std::list<ContractionOperation> ordering;
  std::vector<std::vector<int>> index = {{1, 2}, {3, 4}};
  std::vector<int> cut_values = {5, 6};
  ordering.emplace_back(CutIndex(index, cut_values));
  ASSERT_EQ(ordering.back().op_type, ContractionOperation::CUT);
  EXPECT_EQ(ordering.back().cut.tensors, index);
  EXPECT_EQ(ordering.back().cut.values, cut_values);

  std::vector<int> expand_tensor = {7, 8};
  ordering.emplace_back(ExpandPatch("a", expand_tensor));
  ASSERT_EQ(ordering.back().op_type, ContractionOperation::EXPAND);
  EXPECT_EQ(ordering.back().expand.id, "a");
  EXPECT_EQ(ordering.back().expand.tensor, expand_tensor);

  ordering.emplace_back(MergePatches("a", "b"));
  ASSERT_EQ(ordering.back().op_type, ContractionOperation::MERGE);
  EXPECT_EQ(ordering.back().merge.source_id, "a");
  EXPECT_EQ(ordering.back().merge.target_id, "b");
}

TEST(ContractionTest, RepeatedOperations) {
  // Cannot contract the same tensor twice.
  std::list<ContractionOperation> ordering;
  ordering.emplace_back(ExpandPatch("a", {1, 2}));
  ordering.emplace_back(ExpandPatch("a", {1, 2}));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Cannot cut the same index twice.
  ordering.clear();
  std::vector<std::vector<int>> index = {{1, 2}, {1, 3}};
  ordering.emplace_back(CutIndex(index, {0, 1}));
  ordering.emplace_back(CutIndex(index, {2, 3}));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Cannot merge the same patch twice.
  ordering.clear();
  ordering.emplace_back(MergePatches("a", "b"));
  ordering.emplace_back(MergePatches("a", "c"));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Can merge into the same patch twice.
  ordering.clear();
  ordering.emplace_back(MergePatches("a", "c"));
  ordering.emplace_back(MergePatches("b", "c"));
  EXPECT_TRUE(IsOrderingValid(ordering));
}

TEST(ContractionTest, MergeSafety) {
  // Cannot expand a patch after merging it.
  std::list<ContractionOperation> ordering;
  ordering.emplace_back(MergePatches("a", "b"));
  ordering.emplace_back(ExpandPatch("a", {1, 2}));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Can expand a patch that has been merged into.
  ordering.clear();
  ordering.emplace_back(MergePatches("a", "b"));
  ordering.emplace_back(ExpandPatch("b", {1, 2}));
  EXPECT_TRUE(IsOrderingValid(ordering));
}

TEST(ContractionTest, CutSafety) {
  // Cannot expand the same patch before and after a cut.
  std::list<ContractionOperation> ordering;
  std::vector<std::vector<int>> index = {{1, 2}, {1, 3}};
  std::vector<int> cut_values = {0, 1};
  ordering.emplace_back(ExpandPatch("a", {1, 2}));
  ordering.emplace_back(CutIndex(index, cut_values));
  ordering.emplace_back(ExpandPatch("a", {2, 2}));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Cannot merge into a patch after a cut if it was previously expanded.
  ordering.clear();
  ordering.emplace_back(ExpandPatch("b", {1, 2}));
  ordering.emplace_back(CutIndex(index, cut_values));
  ordering.emplace_back(MergePatches("a", "b"));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Cannot expand a patch after a cut if it was previously merged into.
  ordering.clear();
  ordering.emplace_back(MergePatches("a", "b"));
  ordering.emplace_back(CutIndex(index, cut_values));
  ordering.emplace_back(ExpandPatch("b", {1, 2}));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Cannot merge into a patch before and after a cut.
  ordering.clear();
  ordering.emplace_back(MergePatches("a", "c"));
  ordering.emplace_back(CutIndex(index, cut_values));
  ordering.emplace_back(MergePatches("b", "c"));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Can merge a patch into an empty patch for later use.
  ordering.clear();
  ordering.emplace_back(ExpandPatch("a", {1, 2}));
  ordering.emplace_back(CutIndex(index, cut_values));
  ordering.emplace_back(MergePatches("a", "b"));
  ordering.emplace_back(ExpandPatch("b", {1, 3}));
  EXPECT_TRUE(IsOrderingValid(ordering));
}

// Trivial grid with one "off" qubit and one cut:
//   0 - 1
//   |   x --> cut between (0,1) and (1,1)
//   2 - 3
//       |
//   4   5 --> qubit at (2,0) is off; (2,1) is in final region.
// This circuit should return the input string with amplitude ~= 1 when summing
// over the cut values, but only when the output of (2,1) is a zero.
TEST(ContractionTest, SimpleInitializeData) {
  std::vector<std::vector<Tensor>> tensor_grid;
  for (int i = 0; i < 3; ++i) {
    tensor_grid.push_back(std::vector<Tensor>(2));
  }
  // clang-format off
  std::vector<std::complex<float>> I_4 =
      {1, 0, 0, 0,
       0, 1, 0, 0,
       0, 0, 1, 0,
       0, 0, 0, 1};
  std::vector<std::complex<float>> I_8 =
      {1, 0, 0, 0, 0, 0, 0, 0,
       0, 1, 0, 0, 0, 0, 0, 0,
       0, 0, 1, 0, 0, 0, 0, 0,
       0, 0, 0, 1, 0, 0, 0, 0,
       0, 0, 0, 0, 1, 0, 0, 0,
       0, 0, 0, 0, 0, 1, 0, 0,
       0, 0, 0, 0, 0, 0, 1, 0,
       0, 0, 0, 0, 0, 0, 0, 1};
  std::vector<std::complex<float>> I_2x4 =
      {1, 0, 0, 0,
       0, 1, 0, 0};
  // clang-format on

  tensor_grid[0][0] = Tensor({"(0,0),(0,1)", "(0,0),(1,0)"}, {4, 4}, I_4);
  tensor_grid[0][1] = Tensor({"(0,1),(1,1)", "(0,0),(0,1)"}, {4, 4}, I_4);
  tensor_grid[1][0] = Tensor({"(0,0),(1,0)", "(1,0),(1,1)"}, {4, 4}, I_4);
  tensor_grid[1][1] =
      Tensor({"(0,1),(1,1)", "(1,1),(2,1)", "(1,0),(1,1)"}, {4, 4, 4}, I_8);
  tensor_grid[2][1] = Tensor({"(2,1),(o)", "(1,1),(2,1)"}, {2, 4}, I_2x4);

  std::list<ContractionOperation> ordering;
  ordering.emplace_back(CutIndex({{0, 1}, {1, 1}}));
  ordering.emplace_back(ExpandPatch("a", {0, 1}));
  ordering.emplace_back(ExpandPatch("a", {0, 0}));
  ordering.emplace_back(ExpandPatch("a", {1, 0}));
  ordering.emplace_back(CutIndex({{2, 1}}));
  ordering.emplace_back(ExpandPatch("b", {2, 1}));
  ordering.emplace_back(ExpandPatch("b", {1, 1}));
  ordering.emplace_back(MergePatches("a", "b"));

  std::vector<std::complex<double>> amplitudes(2);
  auto data = ContractionData::Initialize(ordering, &tensor_grid, &amplitudes);

  std::unordered_map<std::string, bool> active_patches;
  for (const auto& patch : data.scratch_list()) {
    active_patches[patch] = false;
  }
  data.ContractGrid(ordering, /*output_index=*/0, active_patches);
  ASSERT_EQ(amplitudes.size(), 2);
  // amplitudes[0] represents <00000|U|00000>, and should return 1.
  EXPECT_FLOAT_EQ(amplitudes[0].real(), 1.0);
  EXPECT_FLOAT_EQ(amplitudes[0].imag(), 0.0);
  // amplitudes[1] represents <00000|U|00001>, and should return 0.
  EXPECT_FLOAT_EQ(amplitudes[1].real(), 0.0);
  EXPECT_FLOAT_EQ(amplitudes[1].imag(), 0.0);
}

// This test demonstrates the creation and validation of a complex contraction
// ordering - specifically, the "alternative" contraction for the 7x7 grid
// defined in this paper: https://arxiv.org/pdf/1811.09599.pdf
TEST(ContractionTest, ExampleOrdering) {
  std::list<ContractionOperation> ordering;
  const std::vector<std::vector<int>> order_A = {
      {0, 0}, {0, 1}, {1, 0}, {1, 1}, {0, 2}, {2, 0}, {1, 2}, {2, 1}, {2, 2}};
  for (const auto& coord : order_A) {
    ordering.emplace_back(ExpandPatch("A", coord));
  }
  const std::vector<std::vector<int>> order_pB = {
      {6, 0}, {5, 0}, {4, 0}, {3, 0}, {6, 1}, {5, 1}, {4, 1}, {3, 1}};
  for (const auto& coord : order_pB) {
    ordering.emplace_back(ExpandPatch("pB", coord));
  }
  const std::vector<std::vector<int>> order_ppD = {
      {6, 6}, {6, 5}, {5, 6}, {5, 5}, {6, 4}, {4, 6}, {5, 4}, {4, 5}, {4, 4}};
  for (const auto& coord : order_ppD) {
    ordering.emplace_back(ExpandPatch("ppD", coord));
  }
  const std::vector<std::vector<std::vector<int>>> cuts_1 = {{{6, 2}, {6, 3}}};
  for (const auto& cut : cuts_1) {
    ordering.emplace_back(CutIndex(cut));
  }
  // Copies tensor "pB" to "B" for reuse.
  ordering.emplace_back(MergePatches("pB", "B"));
  const std::vector<std::vector<int>> order_B = {
      {6, 2}, {5, 2}, {4, 2}, {3, 2}};
  for (const auto& coord : order_B) {
    ordering.emplace_back(ExpandPatch("B", coord));
  }
  ordering.emplace_back(MergePatches("A", "B"));
  // Copies tensor "ppD" to "pD" for reuse.
  ordering.emplace_back(MergePatches("ppD", "pD"));
  const std::vector<std::vector<int>> order_pD = {{6, 3}, {5, 3}, {4, 3}};
  for (const auto& coord : order_pD) {
    ordering.emplace_back(ExpandPatch("pD", coord));
  }
  const std::vector<std::vector<int>> order_pC = {{0, 6}, {1, 6}, {2, 6},
                                                  {0, 5}, {1, 5}, {2, 5}};
  for (const auto& coord : order_pC) {
    ordering.emplace_back(ExpandPatch("pC", coord));
  }
  const std::vector<std::vector<std::vector<int>>> cuts_2 = {{{2, 6}, {3, 6}}};
  for (const auto& cut : cuts_2) {
    ordering.emplace_back(CutIndex(cut));
  }
  // Copies tensor "pD" to "D" for reuse.
  ordering.emplace_back(MergePatches("pD", "D"));
  const std::vector<std::vector<int>> order_D = {
      {3, 6}, {3, 5}, {3, 4}, {3, 3}};
  for (const auto& coord : order_D) {
    ordering.emplace_back(ExpandPatch("D", coord));
  }
  ordering.emplace_back(MergePatches("B", "D"));
  const std::vector<std::vector<int>> order_C = {{0, 4}, {1, 4}, {2, 4},
                                                 {0, 3}, {1, 3}, {2, 3}};
  // These are "terminal cuts" for fast sampling over output values of C.
  for (const auto& tensor : order_C) {
    ordering.emplace_back(CutIndex({tensor}));
  }
  ordering.emplace_back(MergePatches("pC", "C"));
  for (const auto& coord : order_C) {
    ordering.emplace_back(ExpandPatch("C", coord));
  }
  ordering.emplace_back(MergePatches("D", "C"));

  ASSERT_TRUE(IsOrderingValid(ordering));
}

constexpr char kSimpleOrdering[] = R"(# test comment
cut (1,2) 1 3
expand a 1
expand a 0
expand a 2
cut () 5
expand b 5
expand b 3
merge a b
)";
TEST(OrderingParserTest, ParseSimpleOrdering) {
  auto ordering_data = std::stringstream(kSimpleOrdering);
  std::list<ContractionOperation> ordering;
  std::vector<std::vector<int>> qubits_off = {{2, 0}};
  int I = 3;
  int J = 2;
  ASSERT_TRUE(ordering_data_to_contraction_ordering(&ordering_data, I, J,
                                                    qubits_off, &ordering));

  std::list<ContractionOperation> expected_ordering;
  expected_ordering.emplace_back(CutIndex({{0, 1}, {1, 1}}, {1, 2}));
  expected_ordering.emplace_back(ExpandPatch("a", {0, 1}));
  expected_ordering.emplace_back(ExpandPatch("a", {0, 0}));
  expected_ordering.emplace_back(ExpandPatch("a", {1, 0}));
  expected_ordering.emplace_back(CutIndex({{2, 1}}));
  expected_ordering.emplace_back(ExpandPatch("b", {2, 1}));
  expected_ordering.emplace_back(ExpandPatch("b", {1, 1}));
  expected_ordering.emplace_back(MergePatches("a", "b"));

  int op_count = expected_ordering.size();
  ASSERT_EQ(ordering.size(), op_count);
  for (int i = 0; i < op_count; ++i) {
    const auto op = ordering.front();
    const auto expected_op = expected_ordering.front();
    ASSERT_EQ(op.op_type, expected_op.op_type);
    if (op.op_type == ContractionOperation::EXPAND) {
      EXPECT_EQ(op.expand.id, expected_op.expand.id);
      EXPECT_EQ(op.expand.tensor, expected_op.expand.tensor);
    } else if (op.op_type == ContractionOperation::CUT) {
      EXPECT_EQ(op.cut.tensors, expected_op.cut.tensors);
      EXPECT_EQ(op.cut.values, expected_op.cut.values);

    } else if (op.op_type == ContractionOperation::MERGE) {
      EXPECT_EQ(op.merge.source_id, expected_op.merge.source_id);
      EXPECT_EQ(op.merge.target_id, expected_op.merge.target_id);
    }
    ordering.pop_front();
    expected_ordering.pop_front();
  }
}

constexpr char kInvertedCutOrdering[] = R"(# test comment
cut (1,2) 1 0
)";

TEST(OrderingParserTest, ParseCutReordering) {
  auto ordering_data = std::stringstream(kInvertedCutOrdering);
  std::list<ContractionOperation> ordering;
  std::vector<std::vector<int>> qubits_off = {{2, 0}};
  int I = 1;
  int J = 2;
  ASSERT_TRUE(ordering_data_to_contraction_ordering(&ordering_data, I, J,
                                                    qubits_off, &ordering));

  ContractionOperation expected_op(CutIndex({{0, 0}, {0, 1}}, {1, 2}));

  ASSERT_EQ(ordering.size(), 1);
  const auto& op = ordering.front();
  EXPECT_EQ(op.cut.tensors, expected_op.cut.tensors);
  EXPECT_EQ(op.cut.values, expected_op.cut.values);
}

TEST(OrderingParserTest, ParserFailures) {
  std::list<ContractionOperation> ordering;
  std::vector<std::vector<int>> qubits_off = {{2, 0}};
  std::stringstream ordering_data;
  int I = 3;
  int J = 2;

  // Invalid operations cause failures.
  ordering_data = std::stringstream("bad_op 1 2");
  EXPECT_FALSE(ordering_data_to_contraction_ordering(&ordering_data, I, J,
                                                     qubits_off, &ordering));

  // Qubit indices must be within the grid (3x2).
  ordering_data = std::stringstream("expand a 8");
  EXPECT_FALSE(ordering_data_to_contraction_ordering(&ordering_data, I, J,
                                                     qubits_off, &ordering));
  ordering_data = std::stringstream("expand a -2");
  EXPECT_FALSE(ordering_data_to_contraction_ordering(&ordering_data, I, J,
                                                     qubits_off, &ordering));
  ordering_data = std::stringstream("cut () 1 7");
  EXPECT_FALSE(ordering_data_to_contraction_ordering(&ordering_data, I, J,
                                                     qubits_off, &ordering));
  ordering_data = std::stringstream("cut () -1 4");
  EXPECT_FALSE(ordering_data_to_contraction_ordering(&ordering_data, I, J,
                                                     qubits_off, &ordering));

  // Cuts must receive a valid value list.
  ordering_data = std::stringstream("cut 2 3");
  EXPECT_FALSE(ordering_data_to_contraction_ordering(&ordering_data, I, J,
                                                     qubits_off, &ordering));
  // Spaces are not allowed in the value list.
  ordering_data = std::stringstream("cut (1, 2) 2 3");
  EXPECT_FALSE(ordering_data_to_contraction_ordering(&ordering_data, I, J,
                                                     qubits_off, &ordering));
}

TEST(OrderingParserDeathTest, InvalidInput) {
  auto ordering_data = std::stringstream(kSimpleOrdering);
  std::list<ContractionOperation> ordering;
  std::vector<std::vector<int>> qubits_off = {{2, 0}};
  int I = 3;
  int J = 2;

  // Ordering data cannot be null pointer.
  EXPECT_DEATH(ordering_data_to_contraction_ordering(nullptr, I, J, qubits_off,
                                                     &ordering),
               "");

  // Ordering cannot be null pointer.
  EXPECT_DEATH(ordering_data_to_contraction_ordering(&ordering_data, I, J,
                                                     qubits_off, nullptr),
               "");
}

constexpr char kInvalidOrdering[] = R"(# test comment
expand a 1
expand a 1
)";
TEST(OrderingParserDeathTest, InvalidOrderingGenerated) {
  auto ordering_data = std::stringstream(kInvalidOrdering);
  std::list<ContractionOperation> ordering;
  std::vector<std::vector<int>> qubits_off = {{2, 0}};
  int I = 3;
  int J = 2;
  EXPECT_DEATH(ordering_data_to_contraction_ordering(&ordering_data, I, J,
                                                     qubits_off, &ordering),
               "");
}

TEST(ContractionDeathTest, ContractGridInvalidInput) {
  std::list<ContractionOperation> ordering;
  std::vector<std::vector<Tensor>> tensor_grid;
  std::vector<std::complex<double>> amplitudes;

  // Tensor grid cannot be null pointer.
  EXPECT_DEATH(ContractGrid(ordering, nullptr, &amplitudes), "");

  // Amplitudes cannot be null pointer.
  EXPECT_DEATH(ContractGrid(ordering, &tensor_grid, nullptr), "");
}

TEST(ContractionDeathTest, InitializeInvalidInput) {
  std::list<ContractionOperation> ordering;
  std::vector<std::vector<Tensor>> tensor_grid;
  std::vector<std::complex<double>> amplitudes;

  // Tensor grid cannot be null pointer.
  EXPECT_DEATH(ContractionData::Initialize(ordering, nullptr, &amplitudes), "");

  // Amplitudes cannot be null pointer.
  EXPECT_DEATH(ContractionData::Initialize(ordering, &tensor_grid, nullptr),
               "");
}

}  // namespace
}  // namespace qflex

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
