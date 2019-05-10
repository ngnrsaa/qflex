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
  vector<vector<vector<MKLTensor>>> grid_of_tensors;
  s_type scratch[256];

  auto circuit_data = std::stringstream(kNullCircuit);
  circuit_data_to_grid_of_tensors(&circuit_data, 2, 1, 1, "00", "01", {}, {},
                                  grid_of_tensors, scratch);
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

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
