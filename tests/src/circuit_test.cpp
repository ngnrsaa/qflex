#include "circuit.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace qflex {
namespace {

// Passing in an invalid filename to circuit.load().
TEST(CircuitExceptionTest, InvalidFilenameInput) {
  QflexCircuit circuit;
  std::string invalid_filename = "invalid.txt";

  try {
    circuit.load(invalid_filename);
  } catch (std::string msg) {
    EXPECT_THAT(msg,
                testing::HasSubstr("Cannot open circuit file invalid.txt."));
  }
}

TEST(CircuitExceptionTest, PrintGate) {
  QflexGate gate;
  std::stringstream gate_as_string;

  gate.name = "cx";
  gate.cycle = 1;
  gate.qubits = {2, 4};
  gate.params = {3, 5};

  gate_as_string << gate << std::endl;
  EXPECT_EQ(gate_as_string.str(),
            "gate_name: cx\nqubits: 2 4 \nparams: 3 5 \n\n");
}

constexpr char kBadCircuit1[] = R"(0 h 0
0 h 1
9 h 0
9 h 1)";

constexpr char kBadCircuit2[] = R"(6
0 h 0
0 h
9 h 0
9 h 1)";

constexpr char kBadCircuit3[] = R"(2
0 h 0
0 h 1
cz 0 2
9 h 1)";

constexpr char kBadCircuit4[] = R"(2
0 h 0
0 h 1
9 cz 0 2
7 cz 1 3
9 h 1)";

constexpr char kBadCircuit5[] = R"(4
0 h 0
0 h 1
1 t 0
1 t 2
2 fsim(0.1, a) 0 1
4 h 0
4 h 1
)";

constexpr char kBadCircuit6[] = R"(2
0 h 0
0 h 1
1 t 0
1 t a
4 h 0
4 h 1)";

constexpr char kBadCircuit7[] = R"(4
0 h 0
0 h 1
4 t 2
4 t 4 
4 cz 2 4
9 cz 0 2
9 cz 1 3
9 h 1)";

// Testing bad circuits
TEST(CircuitExceptionTest, BadCircuits) {
  QflexCircuit circuit;

  // First line isn't the number of active qubits
  try {
    circuit.load(std::stringstream(kBadCircuit1));
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("[1: 0 h 0] First line in circuit must "
                                        "be the number of active qubits."));
  }

  // Gate is missing parameters.
  try {
    circuit.load(std::stringstream(kBadCircuit2));
  } catch (std::string msg) {
    EXPECT_THAT(msg,
                testing::HasSubstr("[3: 0 h] Gate must be specified as: cycle "
                                   "gate_name[(p1[,p2,...])] q1 [q2, ...]"));
  }

  // First number isn't cycle.
  try {
    circuit.load(std::stringstream(kBadCircuit3));
  } catch (std::string msg) {
    EXPECT_THAT(msg,
                testing::HasSubstr(
                    "[4: cz 0 2] First token must be a valid cycle number."));
  }

  // Cycle isn't increasing.
  try {
    circuit.load(std::stringstream(kBadCircuit4));
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr(
                         "[5: 7 cz 1 3] Cycle number can only increase."));
  }

  // Params aren't numbers.
  try {
    circuit.load(std::stringstream(kBadCircuit5));
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("Params must be valid numbers."));
  }

  // Qubits aren't valid numbers.
  try {
    circuit.load(std::stringstream(kBadCircuit6));
  } catch (std::string msg) {
    EXPECT_THAT(msg,
                testing::HasSubstr("[5: 1 t a] Qubit must be a valid number."));
  }
  // Qubits are being reused.
  try {
    circuit.load(std::stringstream(kBadCircuit7));
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr(
                         "[6: 4 cz 2 4] Qubits can only used once per cycle."));
  }
}

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
14 h 0
14 h 1
14 h 2
14 h 3
14 h 5)";

// Testing loading a simple circuit.
TEST(CircuitTest, SimpleLoadTest) {
  QflexCircuit circuit;
  circuit.load(std::stringstream(kSimpleCircuit));

  EXPECT_EQ(circuit.num_active_qubits, 5);
  EXPECT_EQ(circuit.gates.size(), 22);
}

constexpr char kClearCircuit[] = R"(2
0 h 0
0 h 1
9 h 0
9 h 1)";

// Testing circuit.clear()
TEST(CircuitTest, ClearCircuitTest) {
  QflexCircuit circuit;
  circuit.load(std::stringstream(kClearCircuit));
  circuit.clear();

  EXPECT_EQ(circuit.gates.size(), 0);
  EXPECT_EQ(circuit.num_active_qubits, 0);
}

}  // namespace
}  // namespace qflex

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}