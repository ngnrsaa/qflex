#include "../evaluate_circuit.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace qflex {
namespace {

class GetOutputStatesTest : public testing::Test {
 public:
  void TestOutputExpectations() {
    get_output_states(ordering_, &final_qubits_, &output_states_);
    EXPECT_EQ(final_qubits_, expected_final_qubits_);
    EXPECT_EQ(output_states_, expected_output_states_);
  }

 protected:
  std::vector<std::vector<int>> final_qubits_, expected_final_qubits_;
  std::vector<std::string> output_states_, expected_output_states_;
  std::list<ContractionOperation> ordering_;
};

// ExpandPatch should not affect output states.
TEST_F(GetOutputStatesTest, IgnoresExpandPatch) {
  ordering_.emplace_back(ExpandPatch("a", {0, 0}));
  ordering_.emplace_back(ExpandPatch("a", {0, 1}));
  expected_final_qubits_ = {};
  expected_output_states_ = {""};
  TestOutputExpectations();
}

// MergePatches should not affect output states.
TEST_F(GetOutputStatesTest, IgnoresMergePatches) {
  ordering_.emplace_back(MergePatches("a", "b"));
  ordering_.emplace_back(MergePatches("b", "c"));
  expected_final_qubits_ = {};
  expected_output_states_ = {""};
  TestOutputExpectations();
}

// Non-terminal cuts should not affect output states.
TEST_F(GetOutputStatesTest, IgnoresNonTerminalCuts) {
  ordering_.emplace_back(CutIndex({{0, 0}, {0, 1}}, {1, 2}));
  ordering_.emplace_back(CutIndex({{0, 0}, {1, 0}}, {3}));
  expected_final_qubits_ = {};
  expected_output_states_ = {""};
  TestOutputExpectations();
}

// Terminal cuts are listed in the order applied.
TEST_F(GetOutputStatesTest, TerminalCutsDefineOutputStates) {
  ordering_.emplace_back(CutIndex({{0, 1}}, {0}));
  ordering_.emplace_back(CutIndex({{0, 0}}, {0, 1}));
  ordering_.emplace_back(CutIndex({{1, 0}}, {1}));
  expected_final_qubits_ = {{0, 1}, {0, 0}, {1, 0}};
  expected_output_states_ = {"001", "011"};
  TestOutputExpectations();
}

// Terminal cuts with no values will be evaluated as "0" and "1", since output
// states can only be one of those two values.
TEST_F(GetOutputStatesTest, BlankCutValuesEvaluateBothStates) {
  ordering_.emplace_back(CutIndex({{0, 1}}));
  ordering_.emplace_back(CutIndex({{0, 0}}));
  ordering_.emplace_back(CutIndex({{1, 0}}));
  expected_final_qubits_ = {{0, 1}, {0, 0}, {1, 0}};
  expected_output_states_ = {"000", "001", "010", "011",
                             "100", "101", "110", "111"};
  TestOutputExpectations();
}

// When a mixture of operations are applied, only terminal cuts affect the
// output states.
TEST_F(GetOutputStatesTest, OnlyUseTerminalCuts) {
  ordering_.emplace_back(CutIndex({{0, 1}, {1, 1}}, {1, 2}));
  ordering_.emplace_back(ExpandPatch("a", {0, 1}));
  ordering_.emplace_back(ExpandPatch("a", {0, 0}));
  ordering_.emplace_back(ExpandPatch("a", {1, 0}));
  ordering_.emplace_back(CutIndex({{2, 1}}));
  ordering_.emplace_back(ExpandPatch("b", {2, 1}));
  ordering_.emplace_back(ExpandPatch("b", {1, 1}));
  ordering_.emplace_back(MergePatches("a", "b"));
  expected_final_qubits_ = {{2, 1}};
  expected_output_states_ = {"0", "1"};
  TestOutputExpectations();
}

// Nullptr input in get_output_states()
TEST(GetOutputStatesDeathTest, InvalidInput) {
  std::list<ContractionOperation> ordering;
  std::vector<std::vector<int>> final_qubits;
  std::vector<std::string> output_states;

  // Final qubits cannot be null pointer.
  EXPECT_DEATH(get_output_states(ordering, nullptr, &output_states), "");

  // Output states cannot be null pointer.
  EXPECT_DEATH(get_output_states(ordering, &final_qubits, nullptr), "");
}

// Grid layout with trailing whitespace.
constexpr char kTestGrid[] = R"(0 1 1 0
                                1 1 1 1
                                0 1 0 0
                                )";

TEST(ReadGridTest, ValidGrid3x4) {
  std::stringstream stream(kTestGrid);
  std::vector<std::vector<int>> off_qubits =
      read_grid_layout_from_stream(&stream, 3, 4);
  std::vector<std::vector<int>> expected_off = {
      {0, 0}, {0, 3}, {2, 0}, {2, 2}, {2, 3}};
  EXPECT_EQ(off_qubits, expected_off);
}

TEST(ReadGridTest, ValidGrid6x2) {
  std::stringstream stream(kTestGrid);
  std::vector<std::vector<int>> off_qubits =
      read_grid_layout_from_stream(&stream, 6, 2);
  std::vector<std::vector<int>> expected_off = {
      {0, 0}, {1, 1}, {4, 0}, {5, 0}, {5, 1}};
  EXPECT_EQ(off_qubits, expected_off);
}

// Grid data is too large: 3 * 4 > 5 * 2
TEST(ReadGridDeathTest, InvalidGrid5x2) {
  std::stringstream stream(kTestGrid);
  EXPECT_DEATH(read_grid_layout_from_stream(&stream, 5, 2), "");
}

// Grid data is too small: 3 * 4 < 5 * 3
TEST(ReadGridDeathTest, InvalidGrid5x3) {
  std::stringstream stream(kTestGrid);
  EXPECT_DEATH(read_grid_layout_from_stream(&stream, 5, 3), "");
}

// Below are config strings for a simple grid with one "off" qubit and one cut:
//   0 - 1
//   |   x --> cut between (0,1) and (1,1)
//   2 - 3
//       |
//   4   5 --> qubit at (2,0) is off; (2,1) is in final region.
// This circuit should return the input string with amplitude ~= 1 when summing
// over the cut values, but only when the output of (2,1) is a zero.
constexpr char kSimpleCircuit[] = R"(5
0 h 0
0 h 1
0 h 2
0 h 3
0 h 5
1 cz 0 1
2 cz 0 2
3 cz 1 3
4 cz 2 3
5 cz 3 5
11 cz 0 1
12 cz 0 2
13 cz 1 3
14 cz 2 3
15 cz 3 5
17 h 0
17 h 1
17 h 2
17 h 3
17 h 5)";

constexpr char kSimpleOrdering[] = R"(#
cut () 1 3
expand a 1
expand a 0
expand a 2
cut () 5
expand b 5
expand b 3
merge a b
)";

constexpr char kSimpleGrid[] = R"(1 1
                                  1 1
                                  0 1)";

// Perform a full evaluation of a very simple circuit.
TEST(EvaluateCircuitTest, SimpleCircuit) {
  std::stringstream circuit_data(kSimpleCircuit);
  std::stringstream ordering_data(kSimpleOrdering);
  std::stringstream grid_data(kSimpleGrid);

  QflexInput input;
  input.I = 3;
  input.J = 2;
  input.K = 2;
  input.circuit_data = &circuit_data;
  input.ordering_data = &ordering_data;
  input.grid_data = &grid_data;
  input.initial_state = "00000";
  input.final_state_A = "0000";

  std::vector<std::pair<std::string, std::complex<double>>> amplitudes =
      EvaluateCircuit(&input);
  ASSERT_EQ(amplitudes.size(), 2);
  EXPECT_EQ(amplitudes[0].first, "0000 0");
  EXPECT_NEAR(amplitudes[0].second.real(), 1.0, 1e-5);
  EXPECT_EQ(amplitudes[1].first, "0000 1");
  EXPECT_NEAR(amplitudes[1].second.real(), 0.0, 1e-5);
}

// Nullptr input in EvaluateCircuit()
TEST(EvaluateCircuitDeathTest, InvalidInput) {
  // Input cannot be null pointer.
  EXPECT_DEATH(EvaluateCircuit(nullptr), "");
}

}  // namespace
}  // namespace qflex

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
