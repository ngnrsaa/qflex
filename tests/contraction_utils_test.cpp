#include "../contraction_utils.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

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

TEST(ContractionTest, OperationHandling) {
  ContractionOrdering ordering;
  std::vector<std::vector<int>> index = {{1, 2}, {3, 4}};
  std::vector<int> cut_values = {5, 6};
  ordering.emplace_back(new CutIndex(index, cut_values));
  ASSERT_EQ(ordering.back()->op_type, ContractionOperation::CUT);
  const auto* cut = dynamic_cast<const CutIndex*>(ordering.back().get());
  EXPECT_EQ(cut->tensors, index);
  EXPECT_EQ(cut->values, cut_values);

  std::vector<int> expand_tensor = {7, 8};
  ordering.emplace_back(new ExpandPatch("a", expand_tensor));
  ASSERT_EQ(ordering.back()->op_type, ContractionOperation::EXPAND);
  const auto* expand = dynamic_cast<const ExpandPatch*>(ordering.back().get());
  EXPECT_EQ(expand->id, "a");
  EXPECT_EQ(expand->tensor, expand_tensor);

  ordering.emplace_back(new MergePatches("a", "b"));
  ASSERT_EQ(ordering.back()->op_type, ContractionOperation::MERGE);
  const auto* merge = dynamic_cast<const MergePatches*>(ordering.back().get());
  EXPECT_EQ(merge->source_id, "a");
  EXPECT_EQ(merge->target_id, "b");
}

TEST(ContractionTest, RepeatedOperations) {
  // Cannot contract the same tensor twice.
  ContractionOrdering ordering;
  ordering.emplace_back(new ExpandPatch("a", {1, 2}));
  ordering.emplace_back(new ExpandPatch("a", {1, 2}));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Cannot cut the same index twice.
  ordering.clear();
  std::vector<std::vector<int>> index = {{1, 2}, {1, 3}};
  ordering.emplace_back(new CutIndex(index, {0, 1}));
  ordering.emplace_back(new CutIndex(index, {2, 3}));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Cannot merge the same patch twice.
  ordering.clear();
  ordering.emplace_back(new MergePatches("a", "b"));
  ordering.emplace_back(new MergePatches("a", "c"));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Can merge into the same patch twice.
  ordering.clear();
  ordering.emplace_back(new MergePatches("a", "c"));
  ordering.emplace_back(new MergePatches("b", "c"));
  EXPECT_TRUE(IsOrderingValid(ordering));
}

TEST(ContractionTest, MergeSafety) {
  // Cannot expand a patch after merging it.
  ContractionOrdering ordering;
  ordering.emplace_back(new MergePatches("a", "b"));
  ordering.emplace_back(new ExpandPatch("a", {1, 2}));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Can expand a patch that has been merged into.
  ordering.clear();
  ordering.emplace_back(new MergePatches("a", "b"));
  ordering.emplace_back(new ExpandPatch("b", {1, 2}));
  EXPECT_TRUE(IsOrderingValid(ordering));
}

TEST(ContractionTest, CutSafety) {
  // Cannot expand the same patch before and after a cut.
  ContractionOrdering ordering;
  std::vector<std::vector<int>> index = {{1, 2}, {1, 3}};
  std::vector<int> cut_values = {0, 1};
  ordering.emplace_back(new ExpandPatch("a", {1, 2}));
  ordering.emplace_back(new CutIndex(index, cut_values));
  ordering.emplace_back(new ExpandPatch("a", {2, 2}));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Cannot merge into a patch after a cut if it was previously expanded.
  ordering.clear();
  ordering.emplace_back(new ExpandPatch("b", {1, 2}));
  ordering.emplace_back(new CutIndex(index, cut_values));
  ordering.emplace_back(new MergePatches("a", "b"));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Cannot expand a patch after a cut if it was previously merged into.
  ordering.clear();
  ordering.emplace_back(new MergePatches("a", "b"));
  ordering.emplace_back(new CutIndex(index, cut_values));
  ordering.emplace_back(new ExpandPatch("b", {1, 2}));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Cannot merge into a patch before and after a cut.
  ordering.clear();
  ordering.emplace_back(new MergePatches("a", "c"));
  ordering.emplace_back(new CutIndex(index, cut_values));
  ordering.emplace_back(new MergePatches("b", "c"));
  EXPECT_FALSE(IsOrderingValid(ordering));

  // Can merge a patch into an empty patch for later use.
  ordering.clear();
  ordering.emplace_back(new ExpandPatch("a", {1, 2}));
  ordering.emplace_back(new CutIndex(index, cut_values));
  ordering.emplace_back(new MergePatches("a", "b"));
  ordering.emplace_back(new ExpandPatch("b", {1, 3}));
  EXPECT_TRUE(IsOrderingValid(ordering));
}

// This test demonstrates the creation and validation of a complex contraction
// ordering - specifically, the "alternative" contraction for the 7x7 grid
// defined in this paper: https://arxiv.org/pdf/1811.09599.pdf
TEST(ContractionTest, ExampleOrdering) {
  ContractionOrdering ordering;
  const std::vector<std::vector<int>> order_A = {
      {0, 0}, {0, 1}, {1, 0}, {1, 1}, {0, 2}, {2, 0}, {1, 2}, {2, 1}, {2, 2}};
  for (const auto& coord : order_A) {
    ordering.emplace_back(new ExpandPatch("A", coord));
  }
  const std::vector<std::vector<int>> order_pB = {
      {6, 0}, {5, 0}, {4, 0}, {3, 0}, {6, 1}, {5, 1}, {4, 1}, {3, 1}};
  for (const auto& coord : order_pB) {
    ordering.emplace_back(new ExpandPatch("pB", coord));
  }
  const std::vector<std::vector<int>> order_ppD = {
      {6, 6}, {6, 5}, {5, 6}, {5, 5}, {6, 4}, {4, 6}, {5, 4}, {4, 5}, {4, 4}};
  for (const auto& coord : order_ppD) {
    ordering.emplace_back(new ExpandPatch("ppD", coord));
  }
  const std::vector<std::vector<std::vector<int>>> cuts_1 = {{{6, 2}, {6, 3}}};
  for (const auto& cut : cuts_1) {
    ordering.emplace_back(new CutIndex(cut));
  }
  // Copies tensor "pB" to "B" for reuse.
  ordering.emplace_back(new MergePatches("pB", "B"));
  const std::vector<std::vector<int>> order_B = {
      {6, 2}, {5, 2}, {4, 2}, {3, 2}};
  for (const auto& coord : order_B) {
    ordering.emplace_back(new ExpandPatch("B", coord));
  }
  ordering.emplace_back(new MergePatches("A", "B"));
  // Copies tensor "ppD" to "pD" for reuse.
  ordering.emplace_back(new MergePatches("ppD", "pD"));
  const std::vector<std::vector<int>> order_pD = {{6, 3}, {5, 3}, {4, 3}};
  for (const auto& coord : order_pD) {
    ordering.emplace_back(new ExpandPatch("pD", coord));
  }
  const std::vector<std::vector<int>> order_pC = {{0, 6}, {1, 6}, {2, 6},
                                                  {0, 5}, {1, 5}, {2, 5}};
  for (const auto& coord : order_pC) {
    ordering.emplace_back(new ExpandPatch("pC", coord));
  }
  const std::vector<std::vector<std::vector<int>>> cuts_2 = {{{2, 6}, {3, 6}}};
  for (const auto& cut : cuts_2) {
    ordering.emplace_back(new CutIndex(cut));
  }
  // Copies tensor "pD" to "D" for reuse.
  ordering.emplace_back(new MergePatches("pD", "D"));
  const std::vector<std::vector<int>> order_D = {
      {3, 6}, {3, 5}, {3, 4}, {3, 3}};
  for (const auto& coord : order_D) {
    ordering.emplace_back(new ExpandPatch("D", coord));
  }
  ordering.emplace_back(new MergePatches("B", "D"));
  const std::vector<std::vector<int>> order_C = {{0, 4}, {1, 4}, {2, 4},
                                                 {0, 3}, {1, 3}, {2, 3}};
  // These are "terminal cuts" for fast sampling over output values of C.
  for (const auto& tensor : order_C) {
    ordering.emplace_back(new CutIndex({tensor}));
  }
  ordering.emplace_back(new MergePatches("pC", "C"));
  for (const auto& coord : order_C) {
    ordering.emplace_back(new ExpandPatch("C", coord));
  }
  ordering.emplace_back(new MergePatches("D", "C"));

  ASSERT_TRUE(IsOrderingValid(ordering));
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

}  // namespace

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
