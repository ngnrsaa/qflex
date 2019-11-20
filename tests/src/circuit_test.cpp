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
  gate.name = "cx";
  gate.cycle = 1;
  gate.qubits = {2, 4};
  gate.params = {3, 5};
  // std::cout << gate << std::endl;
}

constexpr char kBadCircuit1[] = R"(
0 h 0
0 h 1
0 h 2)";

TEST(CircuitExceptionTest, BadCircuits) {
  QflexCircuit circuit;
  try {
    circuit.load(std::stringstream(kBadCircuit1));
  // shouldn't this throw an error??
  } catch (std::string msg) {
      std::cout << msg << std::endl;
      EXPECT_THAT(msg, testing::HasSubstr("asdf"));
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

// Testing circuit.clear()
constexpr char kClearCircuit[] = R"(2
0 h 0
0 h 1
9 h 0
9 h 1)";

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