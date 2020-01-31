#include "evaluate_circuit.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace qflex {
namespace {

class GetOutputStatesTest : public testing::Test {
 public:
  void TestOutputExpectations() {
    const auto final_state =
        std::empty(input_.final_states) ? "" : input_.final_states[0];
    const auto final_qubits = get_final_qubits(input_.grid, ordering_);
    const auto output_states = get_output_states(final_state, final_qubits);
    EXPECT_EQ(final_qubits.qubits, expected_final_qubits_);
    EXPECT_EQ(output_states, expected_output_states_);
  }

 protected:
  QflexInput input_;
  std::vector<std::vector<std::size_t>> expected_final_qubits_;
  std::vector<std::string> expected_output_states_;
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

// Terminal cuts are listed in index order, inline with other qubits.
TEST_F(GetOutputStatesTest, TerminalCutsOrderedNormally) {
  input_.grid.I = 3;
  input_.grid.J = 2;
  input_.final_states = {"xx00x0"};
  ordering_.emplace_back(CutIndex({{0, 1}}, {0}));
  ordering_.emplace_back(CutIndex({{0, 0}}, {0, 1}));
  ordering_.emplace_back(CutIndex({{2, 0}}, {1}));
  expected_final_qubits_ = {{0, 1}, {0, 0}, {2, 0}};
  expected_output_states_ = {"000010", "100010"};
  TestOutputExpectations();
}

// Terminal cuts with no values will be evaluated as "0" and "1", since output
// states can only be one of those two values.
TEST_F(GetOutputStatesTest, BlankCutValuesEvaluateBothStates) {
  input_.grid.I = 2;
  input_.grid.J = 2;
  input_.final_states = {"xxx"};
  ordering_.emplace_back(CutIndex({{0, 1}}));
  ordering_.emplace_back(CutIndex({{0, 0}}));
  ordering_.emplace_back(CutIndex({{1, 0}}));
  expected_final_qubits_ = {{0, 1}, {0, 0}, {1, 0}};
  expected_output_states_ = {"000", "001", "100", "101",
                             "010", "011", "110", "111"};
  TestOutputExpectations();
}

// When a mixture of operations are applied, only terminal cuts affect the
// output states.
TEST_F(GetOutputStatesTest, OnlyUseTerminalCuts) {
  input_.grid.I = 3;
  input_.grid.J = 2;
  input_.grid.qubits_off.push_back({2, 0});
  input_.final_states = {"0000x"};
  ordering_.emplace_back(CutIndex({{0, 1}, {1, 1}}, {1, 2}));
  ordering_.emplace_back(ExpandPatch("a", {0, 1}));
  ordering_.emplace_back(ExpandPatch("a", {0, 0}));
  ordering_.emplace_back(ExpandPatch("a", {1, 0}));
  ordering_.emplace_back(CutIndex({{2, 1}}));
  ordering_.emplace_back(ExpandPatch("b", {2, 1}));
  ordering_.emplace_back(ExpandPatch("b", {1, 1}));
  ordering_.emplace_back(MergePatches("a", "b"));
  expected_final_qubits_ = {{2, 1}};
  expected_output_states_ = {"00000", "00001"};
  TestOutputExpectations();
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
1 t 0
1 t 1
1 t 2
1 t 3
1 t 5
2 cz 0 1
3 cx 0 2
4 cx 1 3
5 cz 2 3
6 cz 3 5
11 cz 0 1
12 cx 0 2
13 cx 1 3
14 cz 2 3
15 cx 3 5
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
  input.circuit.load(circuit_data);
  input.ordering.load(ordering_data);
  input.grid.load(grid_data);
  input.initial_states = {"00000"};
  input.final_states = {"1100x"};

  std::vector<std::tuple<std::string, std::string, std::complex<double>>>
      amplitudes = EvaluateCircuit(input);

  ASSERT_EQ(amplitudes.size(), 2ul);
  EXPECT_EQ(std::get<1>(amplitudes[0]), "11000");
  EXPECT_EQ(std::get<1>(amplitudes[1]), "11001");
  EXPECT_NEAR(std::get<2>(amplitudes[0]).real(), 0.10669, 1e-5);
  EXPECT_NEAR(std::get<2>(amplitudes[0]).imag(), 0.04419, 1e-5);
  EXPECT_NEAR(std::get<2>(amplitudes[1]).real(), -0.01831, 1e-5);
  EXPECT_NEAR(std::get<2>(amplitudes[1]).imag(), -0.25758, 1e-5);
}

}  // namespace
}  // namespace qflex

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
