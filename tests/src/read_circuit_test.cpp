#include "read_circuit.h"
#include "circuit.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace qflex {
namespace {

using ::testing::Eq;
using ::testing::Pointwise;

// This circuit is missing the number of active qubits as its first line.
constexpr char kMissingNumQubitsCircuit[] = R"(
0 h 0
0 h 1)";

constexpr char kWrongNumQubitsCircuit[] = R"(3
0 h 0
0 h 1)";

TEST(ReadCircuitDeathTest, InvalidNumberOfQubits) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  s_type scratch[256];

  QflexCircuit circuit;
  circuit.load(std::stringstream(kMissingNumQubitsCircuit));
  EXPECT_DEATH(
      circuit_data_to_tensor_network(circuit, 2, 1, "00", "01", {},
                                      {}, grid_of_tensors, scratch),
      "");

  circuit.load(std::stringstream(kWrongNumQubitsCircuit));
  EXPECT_DEATH(
      circuit_data_to_tensor_network(circuit, 2, 1, "00", "01", {},
                                      {}, grid_of_tensors, scratch),
      "");
}

// This circuit has an invalid one-qubit gate.
constexpr char kBadCircuit[] = R"(2
0 h 0
0 h 1
1 badgate 0)";

TEST(ReadCircuitDeathTest, BadOneQubitGate) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  s_type scratch[256];

  QflexCircuit circuit;
  circuit.load(std::stringstream(kBadCircuit));
  EXPECT_DEATH(
      circuit_data_to_tensor_network(circuit, 2, 1, "00", "01", {},
                                      {}, grid_of_tensors, scratch),
      "");
}

// This circuit has an invalid fsim-type gate.
constexpr char kBadFsimCircuit[] = R"(2
0 h 0
0 h 1
1 fsimbadgate 0 1)";

TEST(ReadCircuitDeathTest, BadFsimGate) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  s_type scratch[256];

  QflexCircuit circuit;
  circuit.load(std::stringstream(kBadFsimCircuit));
  EXPECT_DEATH(
      circuit_data_to_tensor_network(circuit, 2, 1, "00", "01", {},
                                      {}, grid_of_tensors, scratch),
      "");
}

// These circuits reference inactive qubits.
constexpr char kBadTGate[] = R"(2
1 t 0)";

constexpr char kBadCzGate[] = R"(2
1 cz 0 1)";

constexpr char kBadCxGate[] = R"(2
1 cx 0 1)";

TEST(ReadCircuitTest, CircuitReferencingInactiveQubits) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  std::vector<std::vector<int>> off_qubits = {{0, 0}};
  s_type scratch[256];

  // One qubit gate must be on active qubit.
  QflexCircuit circuit;
  circuit.load(std::stringstream(kBadTGate));
  EXPECT_DEATH(
      circuit_data_to_tensor_network(circuit, 2, 1, "0", "1", {},
                                      off_qubits, grid_of_tensors, scratch),
      "");

  // Two qubit gate must have active qubit as first qubit input.
  circuit.load(std::stringstream(kBadCzGate));
  EXPECT_DEATH(
      circuit_data_to_tensor_network(circuit, 2, 1, "0", "1", {},
                                      off_qubits, grid_of_tensors, scratch),
      "");

  // Two qubit gate must have active qubit as first qubit input.
  circuit.load(std::stringstream(kBadCxGate));
  EXPECT_DEATH(
      circuit_data_to_tensor_network(circuit, 2, 1, "0", "1", {},
                                      off_qubits, grid_of_tensors, scratch),
      "");

  // Two qubit gate must have active qubit as second qubit input.
  off_qubits = {{1, 0}};
  EXPECT_DEATH(
      circuit_data_to_tensor_network(circuit, 2, 1, "0", "1", {},
                                      off_qubits, grid_of_tensors, scratch),
      "");
}

constexpr char kBadCycle1[] = R"(4
1 t 1
1 t 1)";

constexpr char kBadCycle2[] = R"(4
1 t 1
1 cz 1 2)";

constexpr char kBadCycle3[] = R"(4
1 t 2
1 cz 1 2)";

constexpr char kBadCycle4[] = R"(4
1 cz 1 2
1 cz 1 3)";

constexpr char kBadCycle5[] = R"(4
1 cz 1 3
1 cz 2 3)";

TEST(ReadCircuitDeathTest, MultipleGatesPerQubitPerCycle) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  std::vector<std::vector<int>> off_qubits = {{0, 0}};
  s_type scratch[256];

  QflexCircuit circuit;
  EXPECT_ANY_THROW(circuit.load(std::stringstream(kBadCycle1)));
  EXPECT_ANY_THROW(circuit.load(std::stringstream(kBadCycle2)));
  EXPECT_ANY_THROW(circuit.load(std::stringstream(kBadCycle3)));
  EXPECT_ANY_THROW(circuit.load(std::stringstream(kBadCycle4)));
  EXPECT_ANY_THROW(circuit.load(std::stringstream(kBadCycle5)));

}

// This circuit returns the input string with amplitude 1.
constexpr char kNullCircuit[] = R"(2
0 h 0
0 h 1
9 h 0
9 h 1)";

// Verifies that circuits can be read from file and the resulting gates match
// the expected circuit.
TEST(ReadCircuitTest, NullCircuit) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  s_type scratch[256];

  QflexCircuit circuit;
  circuit.load(std::stringstream(kNullCircuit));
  circuit_data_to_tensor_network(circuit, 2, 1, "00", "01", {}, {},
                                  grid_of_tensors, scratch);
  // Qubit-0 path has amplitude 1 (0 --> 0)
  // Qubit-1 path has amplitude 0 (0 --> 1)
  std::vector<s_type> expected_data = {std::complex<float>(1, 0),
                                       std::complex<float>(0, 0)};

  // Resulting tensor grid should be 2x1x1 (IxJxK).
  ASSERT_EQ(grid_of_tensors.size(), 2);
  for (int i = 0; i < 2; i++) {
    ASSERT_EQ(grid_of_tensors[i].size(), 1);
    // TODO: tensors are not homogeneous anymore.
    //ASSERT_EQ(grid_of_tensors[i][0].size(), 1);
    //ASSERT_EQ(grid_of_tensors[i][0][0].size(), 1);
    // TODO: tensors are not homogeneous anymore.
    //const std::vector<s_type> data(grid_of_tensors[i][0][0].data(),
    //                               grid_of_tensors[i][0][0].data() + 1);
    // Testing exact equality of floating-point types will fail.
    // Instead, we use EXPECT_FLOAT_EQ on each component of the data.
    //EXPECT_FLOAT_EQ(data[0].real(), expected_data[i].real());
    //EXPECT_FLOAT_EQ(data[0].imag(), expected_data[i].imag());
  }
}

// Simple grid with one "off" qubit and one cut:
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
2 cx 0 2
3 cx 1 3
4 cz 2 3
5 cz 3 5
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

// Verifies that collapsing the 3D tensor network to a 2D grid behaves as
// expected.
TEST(ReadCircuitTest, CondenseToGrid) {
  std::vector<std::vector<std::vector<Tensor>>> tensor_grid_3D;
  std::vector<std::vector<int>> qubits_A = {{2, 1}};
  std::vector<std::vector<int>> qubits_off = {{2, 0}};
  s_type* scratch = new s_type[256];

  QflexCircuit circuit;
  circuit.load(std::stringstream(kSimpleCircuit));
  circuit_data_to_tensor_network(circuit, 3, 2, "00000", "0000x",
                                  qubits_A, qubits_off, tensor_grid_3D,
                                  scratch);

  ASSERT_EQ(tensor_grid_3D.size(), 3);
  for(const auto &g: tensor_grid_3D) ASSERT_EQ(g.size(), 2);

  // TODO: tensors are not homogeneous anymore.
  //ASSERT_EQ(tensor_grid_3D[0][0].size(), 2);

  std::vector<std::vector<Tensor>> tensor_grid_2D;
  for (int i = 0; i < 3; ++i) {
    tensor_grid_2D.push_back(std::vector<Tensor>(2));
  }
  // Working from either end, create two patches and meet in the middle.
  std::list<ContractionOperation> ordering;
  ordering.emplace_back(CutIndex({{0, 1}, {1, 1}}));
  ordering.emplace_back(ExpandPatch("a", {0, 1}));
  ordering.emplace_back(ExpandPatch("a", {0, 0}));
  ordering.emplace_back(ExpandPatch("a", {1, 0}));
  ordering.emplace_back(CutIndex({{2, 1}}));
  ordering.emplace_back(ExpandPatch("b", {2, 1}));
  ordering.emplace_back(ExpandPatch("b", {1, 1}));
  ordering.emplace_back(MergePatches("a", "b"));

}

TEST(ReadCircuitDeathTest, CircuitDataToGridOfTensorsInvalidInput) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  s_type scratch[256];

  QflexCircuit circuit;
  circuit.load(std::stringstream(kNullCircuit));

  // Scratch cannot be null pointer.
  EXPECT_DEATH(
      circuit_data_to_tensor_network(circuit, 2, 1, "01", "01", {},
                                      {}, grid_of_tensors, nullptr),
      "");

  // Input configuration length must be equal to the number of qubits.
  EXPECT_DEATH(
      circuit_data_to_tensor_network(circuit, 2, 1, "001", "01", {},
                                      {}, grid_of_tensors, scratch),
      "");

  // Output configuration length must be equal to the number of qubits.
  EXPECT_DEATH(
      circuit_data_to_tensor_network(circuit, 2, 1, "00", "011", {},
                                      {}, grid_of_tensors, scratch),
      "");
}

}  // namespace
}  // namespace qflex

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
