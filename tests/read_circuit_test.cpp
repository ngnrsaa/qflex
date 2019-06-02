#include "../read_circuit.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {

using ::testing::Eq;
using ::testing::Pointwise;

// This circuit returns the input string with amplitude 1.
constexpr char kNullCircuit[] = R"(2
0 h 0
0 h 1
9 h 0
9 h 1)";

// Verifies that circuits can be read from file and the resulting gates match
// the expected circuit.
TEST(ReadCircuitTest, NullCircuit) {
  std::vector<std::vector<std::vector<MKLTensor>>> grid_of_tensors;
  s_type scratch[256];

  auto circuit_data = std::stringstream(kNullCircuit);
  internal::circuit_data_to_grid_of_tensors(&circuit_data, 2, 1, 1, "00", "01",
                                            {}, {}, grid_of_tensors, scratch);
  // Qubit-0 path has amplitude 1 (0 --> 0)
  // Qubit-1 path has amplitude 0 (0 --> 1)
  std::vector<s_type> expected_data = {std::complex<float>(1, 0),
                                       std::complex<float>(0, 0)};

  // Resulting tensor grid should be 2x1x1 (IxJxK).
  ASSERT_EQ(grid_of_tensors.size(), 2);
  for (int i = 0; i < 2; i++) {
    ASSERT_EQ(grid_of_tensors[i].size(), 1);
    ASSERT_EQ(grid_of_tensors[i][0].size(), 1);
    ASSERT_EQ(grid_of_tensors[i][0][0].size(), 1);
    const std::vector<s_type> data(grid_of_tensors[i][0][0].data(),
                                   grid_of_tensors[i][0][0].data() + 1);
    // Testing exact equality of floating-point types will fail.
    // Instead, we use EXPECT_FLOAT_EQ on each component of the data.
    EXPECT_FLOAT_EQ(data[0].real(), expected_data[i].real());
    EXPECT_FLOAT_EQ(data[0].imag(), expected_data[i].imag());
  }
}

// Simple grid with one 'off' qubit and one cut:
//   0 - 1
//   |   x --> cut between (0,1) and (1,1)
//   2 - 3
//       |
//   4   5 --> qubit at (2,0) is off; (2,1) is in final region.
// This circuit should return the input string with amplitude ~= 1 when summing
// over the cut values.
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

// Verifies that collapsing the 3D tensor network to a 2D grid behaves as
// expected.
TEST(ReadCircuitTest, CondenseToGrid) {
  std::vector<std::vector<std::vector<MKLTensor>>> tensor_grid_3D;
  std::vector<std::vector<int>> qubits_A = {{2, 1}};
  std::vector<std::vector<int>> qubits_off = {{2, 0}};
  s_type * scratch = new s_type[256];
  auto circuit_data = std::stringstream(kSimpleCircuit);
  internal::circuit_data_to_grid_of_tensors(&circuit_data, 3, 2, 2, "00000",
                                            "0000", qubits_A, qubits_off,
                                            tensor_grid_3D, scratch);
  ASSERT_EQ(tensor_grid_3D.size(), 3);
  ASSERT_EQ(tensor_grid_3D[0].size(), 2);
  ASSERT_EQ(tensor_grid_3D[0][0].size(), 2);

  std::vector<std::vector<MKLTensor>> tensor_grid_2D;
  for (int i = 0; i < 3; ++i) {
    tensor_grid_2D.push_back(std::vector<MKLTensor>(2));
  }
  // Working from either end, create two patches and meet in the middle.
  ContractionOrdering ordering;
  ordering.emplace_back(new ExpandPatch('a', {0, 1}));
  ordering.emplace_back(new ExpandPatch('a', {0, 0}));
  ordering.emplace_back(new ExpandPatch('a', {1, 0}));
  ordering.emplace_back(new CutIndex({0, 1}, {1, 1}));
  ordering.emplace_back(new ExpandPatch('b', {2, 1}));
  ordering.emplace_back(new ExpandPatch('b', {1, 1}));
  ordering.emplace_back(new MergePatches('a', 'b'));

  grid_of_tensors_3D_to_2D(tensor_grid_3D, tensor_grid_2D, qubits_A, qubits_off,
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

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
