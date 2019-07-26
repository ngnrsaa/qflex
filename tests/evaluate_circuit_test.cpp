#include "../evaluate_circuit.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace qflex {
namespace {

class GetOutputStatesTest : public testing::Test {
 public:
  bool TestOutputExpectations() {
    get_output_states(ordering_, &final_qubits_, &output_states_);
    EXPECT_EQ(final_qubits_, expected_final_qubits_);
    EXPECT_EQ(output_states_, expected_output_states_);
  }

 protected:
  std::vector<std::vector<int>> final_qubits_, expected_final_qubits_;
  std::vector<std::string> output_states_, expected_output_states_;
  ContractionOrdering ordering_;
};

// ExpandPatch should not affect output states.
TEST_F(GetOutputStatesTest, IgnoresExpandPatch) {
  ordering_.emplace_back(new ExpandPatch("a", {0, 0}));
  ordering_.emplace_back(new ExpandPatch("a", {0, 1}));
  expected_final_qubits_ = {};
  expected_output_states_ = {""};
  TestOutputExpectations();
}

// MergePatches should not affect output states.
TEST_F(GetOutputStatesTest, IgnoresMergePatches) {
  ordering_.emplace_back(new MergePatches("a", "b"));
  ordering_.emplace_back(new MergePatches("b", "c"));
  expected_final_qubits_ = {};
  expected_output_states_ = {""};
  TestOutputExpectations();
}

// Non-terminal cuts should not affect output states.
TEST_F(GetOutputStatesTest, IgnoresNonTerminalCuts) {
  ordering_.emplace_back(new CutIndex({{0, 0}, {0, 1}}, {1, 2}));
  ordering_.emplace_back(new CutIndex({{0, 0}, {1, 0}}, {3}));
  expected_final_qubits_ = {};
  expected_output_states_ = {""};
  TestOutputExpectations();
}

// Terminal cuts are listed in the order applied.
TEST_F(GetOutputStatesTest, TerminalCutsDefineOutputStates) {
  ordering_.emplace_back(new CutIndex({{0, 1}}, {0}));
  ordering_.emplace_back(new CutIndex({{0, 0}}, {0, 1}));
  ordering_.emplace_back(new CutIndex({{1, 0}}, {1}));
  expected_final_qubits_ = {{0, 1}, {0, 0}, {1, 0}};
  expected_output_states_ = {"001", "011"};
  TestOutputExpectations();
}

// Terminal cuts with no values will be evaluated as "0" and "1", since output
// states can only be one of those two values.
TEST_F(GetOutputStatesTest, BlankCutValuesEvaluateBothStates) {
  ordering_.emplace_back(new CutIndex({{0, 1}}));
  ordering_.emplace_back(new CutIndex({{0, 0}}));
  ordering_.emplace_back(new CutIndex({{1, 0}}));
  expected_final_qubits_ = {{0, 1}, {0, 0}, {1, 0}};
  expected_output_states_ = {"000", "001", "010", "011",
                             "100", "101", "110", "111"};
  TestOutputExpectations();
}

// When a mixture of operations are applied, only terminal cuts affect the
// output states.
TEST_F(GetOutputStatesTest, OnlyUseTerminalCuts) {
  ordering_.emplace_back(new CutIndex({{0, 1}, {1, 1}}, {1, 2}));
  ordering_.emplace_back(new ExpandPatch("a", {0, 1}));
  ordering_.emplace_back(new ExpandPatch("a", {0, 0}));
  ordering_.emplace_back(new ExpandPatch("a", {1, 0}));
  ordering_.emplace_back(new CutIndex({{2, 1}}));
  ordering_.emplace_back(new ExpandPatch("b", {2, 1}));
  ordering_.emplace_back(new ExpandPatch("b", {1, 1}));
  ordering_.emplace_back(new MergePatches("a", "b"));
  expected_final_qubits_ = {{2, 1}};
  expected_output_states_ = {"0", "1"};
  TestOutputExpectations();
}

// TODO(martinop): add tests for reading grid and full circuit evaluation.

}  // namespace
}  // namespace qflex

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
