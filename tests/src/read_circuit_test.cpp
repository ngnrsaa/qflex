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

TEST(ReadCircuitExceptionTest, InvalidNumberOfQubits) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  s_type scratch[256];

  QflexCircuit circuit;
  circuit.load(std::stringstream(kMissingNumQubitsCircuit));
  try {
    circuit_data_to_tensor_network(circuit, 2, 1, "00", "01", {}, {},
                                   grid_of_tensors, scratch);
    FAIL()
        << "Expected circuit_data_to_tensor_network() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg, testing::HasSubstr(
                 "The number of active qubits read from the file: 0, does not "
                 "match the number of active qubits read from the grid: 2."));
  }

  circuit.load(std::stringstream(kWrongNumQubitsCircuit));
  try {
    circuit_data_to_tensor_network(circuit, 2, 1, "00", "01", {}, {},
                                   grid_of_tensors, scratch);
    FAIL()
        << "Expected circuit_data_to_tensor_network() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg, testing::HasSubstr(
                 "The number of active qubits read from the file: 3, does not "
                 "match the number of active qubits read from the grid: 2."));
  }
}

// This circuit has an invalid one-qubit gate.
constexpr char kBadCircuit[] = R"(2
0 h 0
0 h 1
1 badgate 0)";

TEST(ReadCircuitExceptionTest, BadOneQubitGate) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  s_type scratch[256];

  QflexCircuit circuit;
  circuit.load(std::stringstream(kBadCircuit));
  try {
    circuit_data_to_tensor_network(circuit, 2, 1, "00", "01", {}, {},
                                   grid_of_tensors, scratch);
    FAIL()
        << "Expected circuit_data_to_tensor_network() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("Invalid gate name provided: badgate"));
  }
}

// This circuit has an invalid fsim-type gate.
constexpr char kBadFsimCircuit[] = R"(2
0 h 0
0 h 1
1 fsimbadgate 0 1)";

TEST(ReadCircuitExceptionTest, BadFsimGate) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  s_type scratch[256];

  QflexCircuit circuit;
  circuit.load(std::stringstream(kBadFsimCircuit));
  try {
    circuit_data_to_tensor_network(circuit, 2, 1, "00", "01", {}, {},
                                   grid_of_tensors, scratch);
    FAIL()
        << "Expected circuit_data_to_tensor_network() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg,
                testing::HasSubstr("Invalid gate name provided: fsimbadgate"));
  }
}

// These circuits reference inactive qubits.
constexpr char kBadTGate[] = R"(1
1 t 0)";

constexpr char kBadCzGate[] = R"(1
1 cz 0 1)";

// Broken test
TEST(ReadCircuitTest, CircuitReferencingInactiveQubits) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  std::vector<std::vector<int>> off_qubits = {{0, 0}};
  s_type scratch[256];

  // One qubit gate must be on active qubit.
  QflexCircuit circuit;
  circuit.load(std::stringstream(kBadTGate));
  try {
    circuit_data_to_tensor_network(circuit, 2, 1, "0", "1", {}, off_qubits,
                                   grid_of_tensors, scratch);
    FAIL()
        << "Expected circuit_data_to_tensor_network() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg, testing::HasSubstr(
                 "Qubit '0' in gate '1 t 0' must correspond to an active qubit."));
  }

  // Two qubit gate must have active qubit as first qubit input.
  circuit.load(std::stringstream(kBadCzGate));
  try {
    circuit_data_to_tensor_network(circuit, 2, 1, "0", "1", {}, off_qubits,
                                   grid_of_tensors, scratch);
    FAIL()
        << "Expected circuit_data_to_tensor_network() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg,
        testing::HasSubstr(
            "Qubit '0' in gate '1 cz 0 1' must correspond to an active qubit."));
  }

  // Two qubit gate must have active qubit as second qubit input.
  off_qubits = {{1, 0}};
  try {
    circuit_data_to_tensor_network(circuit, 2, 1, "0", "1", {}, off_qubits,
                                   grid_of_tensors, scratch);
    FAIL()
        << "Expected circuit_data_to_tensor_network() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg, testing::HasSubstr(
                 "Qubit '1' in gate '1 cz 0 1' must correspond to an active qubit"));
  }

  // Off qubit outside the grid
  off_qubits = {{0, 1}};
  try {
    circuit_data_to_tensor_network(circuit, 2, 1, "0", "1", {}, off_qubits,
                                   grid_of_tensors, scratch);
    FAIL()
        << "Expected circuit_data_to_tensor_network() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg, testing::HasSubstr(
                 "Off qubit '(0, 1)' is outside the grid."));
  }
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

TEST(ReadCircuitExceptionTest, MultipleGatesPerQubitPerCycle) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  std::vector<std::vector<int>> off_qubits = {{0, 0}};
  s_type scratch[256];

  // TODO: fix tests, currently these tests do not test if a qubit is
  // being reused in a cycle.
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
  circuit_data_to_tensor_network(circuit, 2, 1, "00", "10", {}, {},
                                 grid_of_tensors, scratch);
  // Qubit-0 path has amplitude 1 (0 --> 0)
  // Qubit-1 path has amplitude 0 (0 --> 1)
  std::vector<s_type> expected_data = {std::complex<float>(1, 0),
                                       std::complex<float>(0, 0)};

  // Resulting tensor grid should be 2x1 (IxJ) ..
  ASSERT_EQ(std::size(grid_of_tensors), 2);
  for (std::size_t i = 0; i < 2; ++i) {
    const auto &tensor = grid_of_tensors[i];
    ASSERT_EQ(std::size(tensor), 1);

    // .. and each qubit tensor should have 4 tensors ..
    for (const auto &t : tensor) {
      ASSERT_EQ(std::size(t), 4);

      // .. of size 1, 2, 2, 1 (x2, because of complex numbers) respectively
      ASSERT_EQ(std::size(t[0]), 1 * 2);
      ASSERT_EQ(std::size(t[1]), 2 * 2);
      ASSERT_EQ(std::size(t[2]), 2 * 2);
      ASSERT_EQ(std::size(t[3]), 1 * 2);

      // Each qubit tensor should have 4 elements (delta_0, h, h, delta_0)
      // and (delta_0, h, h, delta_1) respectively
      ASSERT_EQ(t[0].data()[0], std::complex<float>(1, 0));
      ASSERT_EQ(t[1].data()[0], std::complex<float>(1. / M_SQRT2, 0));
      ASSERT_EQ(t[1].data()[1], std::complex<float>(1. / M_SQRT2, 0));
      ASSERT_EQ(t[2].data()[0], std::complex<float>(1. / M_SQRT2, 0));
      ASSERT_EQ(t[2].data()[1], std::complex<float>(1. / M_SQRT2, 0));
      ASSERT_EQ(t[3].data()[0], std::complex<float>(i, 0));
    }
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
  s_type *scratch = new s_type[256];

  QflexCircuit circuit;
  circuit.load(std::stringstream(kSimpleCircuit));
  circuit_data_to_tensor_network(circuit, 3, 2, "00000", "0000x", qubits_A,
                                 qubits_off, tensor_grid_3D, scratch);

  ASSERT_EQ(tensor_grid_3D.size(), 3);
  for (const auto &tensor : tensor_grid_3D) ASSERT_EQ(tensor.size(), 2);

  // TODO Add check on tensors

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

  flatten_grid_of_tensors(tensor_grid_3D, tensor_grid_2D, qubits_A, qubits_off,
                          ordering, scratch);

  // Verify that index ordering follows this pattern:
  //   1) Final-region indices ("<index>,(o)")
  //   2) Cut indices
  //   3) Same-patch indices in order by index
  //   4) Cross-patch indices in order by patch
  std::vector<std::vector<std::vector<std::string>>> expected_indices = {
      {
          // qubit (0,0) - normal
          {"(0,0),(0,1)", "(0,0),(1,0)"},
          // qubit (0,1) - cut index
          {"(0,1),(1,1)", "(0,0),(0,1)"},
      },
      {
          // qubit (1,0) - cross-patch index
          {"(0,0),(1,0)", "(1,0),(1,1)"},
          // qubit (1,1) - cut index and cross-patch index
          {"(0,1),(1,1)", "(1,1),(2,1)", "(1,0),(1,1)"},
      },
      {
          // qubit (2,0) - off
          {},
          // qubit (2,1) - final-region index
          {"(2,1),(o)", "(1,1),(2,1)"},
      },
  };
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_EQ(tensor_grid_2D[i][j].get_indices(), expected_indices[i][j]);
    }
  }
}

TEST(ReadCircuitExceptionTest, CircuitDataToGridOfTensorsInvalidInput) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors;
  s_type scratch[256];

  QflexCircuit circuit;
  circuit.load(std::stringstream(kNullCircuit));

  // Scratch cannot be null pointer.
  try {
    circuit_data_to_tensor_network(circuit, 2, 1, "01", "01", {}, {},
                                   grid_of_tensors, nullptr);
    FAIL()
        << "Expected circuit_data_to_tensor_network() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("Scratch must be non-null."));
  }

  // Input configuration length must be equal to the number of qubits.
  try {
    circuit_data_to_tensor_network(circuit, 2, 1, "001", "01", {}, {},
                                   grid_of_tensors, scratch);
    FAIL()
        << "Expected circuit_data_to_tensor_network() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("Size of initial_conf: 3, must be "
                                        "equal to the number of qubits: 2."));
  }

  // Output configuration length must be equal to the number of qubits.
  try {
    circuit_data_to_tensor_network(circuit, 2, 1, "01", "011", {}, {},
                                   grid_of_tensors, scratch);
    FAIL()
        << "Expected circuit_data_to_tensor_network() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("Size of final_conf: 3, must be equal "
                                        "to size of initial_conf: 2."));
  }
}

TEST(ReadCircuitExceptionTest, GridOfTensors3DTo2DInvalidInput) {
  std::vector<std::vector<std::vector<Tensor>>> grid_of_tensors_3D;
  std::vector<std::vector<Tensor>> grid_of_tensors_2D;
  std::optional<std::vector<std::vector<int>>> A;
  std::optional<std::vector<std::vector<int>>> off;
  std::list<ContractionOperation> ordering;

  // Scratch cannot be null pointer.
  try {
    flatten_grid_of_tensors(grid_of_tensors_3D, grid_of_tensors_2D, A, off,
                            ordering, nullptr);
    FAIL() << "Expected flatten_grid_of_tensors() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("Scratch must be non-null."));
  }
}

}  // namespace
}  // namespace qflex

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
