#include "../mkl_tensor.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {

using ::testing::Eq;
using ::testing::Pointwise;

// Creates an empty tensor and runs basic sanity checks on it.
TEST(MKLTensorTest, EmptyTensor) {
  std::vector<std::string> indices = {"a", "b"};
  std::vector<size_t> dimensions = {2, 4};

  // Automatically allocates new space.
  MKLTensor tensor(indices, dimensions);
  ASSERT_EQ(tensor.get_indices(), indices);
  ASSERT_EQ(tensor.get_dimensions(), dimensions);
  ASSERT_EQ(tensor.size(), 8);
  ASSERT_EQ(tensor.num_zeros(), 8);

  // Verify that data is initialized to zero.
  std::vector<std::complex<float>> read_data(tensor.data(),
                                             tensor.data() + tensor.size());
  for (int i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), 0);
    ASSERT_FLOAT_EQ(read_data[i].imag(), 0);
  }

  const auto dict = tensor.get_index_to_dimension();
  ASSERT_EQ(dict.at("a"), 2);
  ASSERT_EQ(dict.at("b"), 4);

  // Test index renaming as well.
  tensor.rename_index("a", "c");
  std::vector<std::string> expected_indices = {"c", "b"};
  ASSERT_EQ(tensor.get_indices(), expected_indices);
}

// Loads a tensor from data and runs basic sanity checks on it.
TEST(MKLTensorTest, LoadData) {
  std::vector<std::string> indices = {"a", "b"};
  std::vector<size_t> dimensions = {2, 2};
  std::vector<std::complex<float>> data = {
      std::complex<float>(0, 1), std::complex<float>(2, 3),
      std::complex<float>(4, 5), std::complex<float>(6, 7),
  };
  MKLTensor tensor(indices, dimensions, data);
  ASSERT_EQ(tensor.get_indices(), indices);
  ASSERT_EQ(tensor.get_dimensions(), dimensions);
  ASSERT_EQ(tensor.size(), 4);
  ASSERT_EQ(tensor.num_zeros(), 0);

  // Verify that data is read in from vector as expected.
  std::vector<std::complex<float>> read_data(tensor.data(),
                                             tensor.data() + tensor.size());
  for (int i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), data[i].real());
    ASSERT_FLOAT_EQ(read_data[i].imag(), data[i].imag());
  }

  // Test scalar multiplication as well.
  tensor.scalar_multiply(std::complex<float>(0, 1));
  read_data = std::vector<std::complex<float>>(tensor.data(),
                                               tensor.data() + tensor.size());
  for (int i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), -1 * data[i].imag());
    ASSERT_FLOAT_EQ(read_data[i].imag(), data[i].real());
  }
}

// Projects a tensor onto a single value of an index and verifies that the
// output tensor only contains data from that slice of the original tensor.
TEST(MKLTensorTest, TensorProjection) {
  std::vector<std::string> indices = {"a", "b", "c"};
  std::vector<size_t> dimensions = {2, 2, 2};
  std::vector<std::complex<float>> data;
  for (int i = 0; i < 8; i++) {
    data.push_back(std::complex<float>(i, 0));
  }

  MKLTensor tensor(indices, dimensions, data);
  std::vector<std::string> expected_indices = {"b", "c"};
  std::vector<size_t> expected_dimensions = {2, 2};

  MKLTensor projection_tensor_1({"x", "y"}, {2, 2});
  tensor.project("a", 1, projection_tensor_1);
  ASSERT_EQ(projection_tensor_1.get_indices(), expected_indices);
  ASSERT_EQ(projection_tensor_1.get_dimensions(), expected_dimensions);

  MKLTensor projection_tensor_2({""}, {64});
  tensor.project("a", 1, projection_tensor_2);
  ASSERT_EQ(projection_tensor_2.get_indices(), expected_indices);
  ASSERT_EQ(projection_tensor_2.get_dimensions(), expected_dimensions);

  const int psize = projection_tensor_1.size();

  const std::vector<std::complex<float>> proj_data_1(
      projection_tensor_1.data(), projection_tensor_1.data() + psize);
  const std::vector<std::complex<float>> proj_data_2(
      projection_tensor_2.data(), projection_tensor_1.data() + psize);
  for (int i = 0; i < psize; ++i) {
    ASSERT_FLOAT_EQ(proj_data_1[i].real(), data[i + psize].real());
    ASSERT_FLOAT_EQ(proj_data_2[i].real(), data[i + psize].real());
  }
}

// Bundles indices of a tensor and verifies that dimensions change while data
// remains unaffected.
TEST(MKLTensorTest, IndexBundling) {
  std::vector<std::string> indices = {"a", "b", "c", "d"};
  std::vector<size_t> dimensions = {2, 2, 2, 2};
  std::vector<std::complex<float>> data;
  for (int i = 0; i < 16; i++) {
    data.push_back(std::complex<float>(i, 0));
  }

  MKLTensor tensor(indices, dimensions, data);
  tensor.bundle({"a", "b", "c"}, "abc");
  std::vector<std::string> expected_indices = {"abc", "d"};
  std::vector<size_t> expected_dimensions = {8, 2};
  ASSERT_EQ(tensor.get_indices(), expected_indices);
  ASSERT_EQ(tensor.get_dimensions(), expected_dimensions);

  std::vector<std::complex<float>> read_data(tensor.data(),
                                             tensor.data() + tensor.size());
  for (int i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), data[i].real());
  }
}

// Reorders indices of a tensor and verifies that data changes accordingly.
TEST(MKLTensorTest, IndexReordering) {
  std::vector<std::string> indices = {"a", "b", "c"};
  std::vector<size_t> dimensions = {2, 2, 2};
  std::vector<std::complex<float>> data;
  for (int i = 0; i < 8; i++) {
    data.push_back(std::complex<float>(i, 0));
  }

  MKLTensor tensor(indices, dimensions, data);
  std::vector<std::string> expected_indices = {"b", "c", "a"};
  std::array<std::complex<float>, 8> scratch;
  tensor.reorder(expected_indices, scratch.data());
  ASSERT_EQ(tensor.get_indices(), expected_indices);
  ASSERT_EQ(tensor.get_dimensions(), dimensions);

  std::vector<std::complex<float>> read_data(tensor.data(),
                                             tensor.data() + tensor.size());
  std::vector<std::complex<float>> expected_data = {
    std::complex<float>(0, 0), std::complex<float>(4, 0),
    std::complex<float>(1, 0), std::complex<float>(5, 0),
    std::complex<float>(2, 0), std::complex<float>(6, 0),
    std::complex<float>(3, 0), std::complex<float>(7, 0),
  };
  for (int i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), expected_data[i].real());
  }
}

// Multiplies two tensors and verify shape, indices, and data of the result.
TEST(MKLTensorTest, Multiply) {
  std::vector<std::string> indices_a = {"a", "b", "c"};
  std::vector<size_t> dimensions_a = {2, 2, 2};
  std::vector<std::complex<float>> data_a;
  for (int i = 0; i < 8; i++) {
    data_a.push_back(std::complex<float>(i, 0));
  }
  MKLTensor tensor_a(indices_a, dimensions_a, data_a);

  std::vector<std::string> indices_b = {"b", "c", "d"};
  std::vector<size_t> dimensions_b = {2, 2, 2};
  std::vector<std::complex<float>> data_b;
  for (int i = 8; i > 0; i--) {
    data_b.push_back(std::complex<float>(i, 0));
  }
  MKLTensor tensor_b(indices_b, dimensions_b, data_b);

  std::vector<std::string> indices_c = {"x"};
  std::vector<size_t> dimensions_c = {16};
  MKLTensor tensor_c(indices_c, dimensions_c);

  std::array<std::complex<float>, 16> scratch;
  multiply(tensor_a, tensor_b, tensor_c, scratch.data());
  std::vector<std::string> expected_indices = {"a", "d"};
  std::vector<size_t> expected_dimensions = {2, 2};
  ASSERT_EQ(tensor_c.get_indices(), expected_indices);
  ASSERT_EQ(tensor_c.get_dimensions(), expected_dimensions);

  std::vector<std::complex<float>> expected_data = {20, 14, 100, 78};
  std::vector<std::complex<float>> read_data(tensor_c.data(),
                                             tensor_c.data() + tensor_c.size());
  for (int i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), expected_data[i].real());
  }
}

// Verifies that a tensor retains its initialized capacity.
TEST(MKLTensorDeathTest, Capacity) {
  std::vector<std::string> indices = {"a", "b"};
  std::vector<size_t> dimensions = {2, 4};

  // Create an uninitialized tensor.
  MKLTensor tensor;

  // Allocate 64 units of space, which should always be available.
  tensor = MKLTensor({""}, {64});

  // Reduce size to 16 units; capacity is still 64 units.
  tensor.set_indices_and_dimensions({"a", "b"}, {4, 4});

  // Increase size to 32 units; capacity is still 64 units.
  tensor.set_indices_and_dimensions({"x"}, {32});

  // Reduce size to 2 units; capacity is still 64 units.
  tensor.set_indices_and_dimensions({"z"}, {2});

  // Increase size to 8 units; capacity is still 64 units.
  tensor.set_indices_and_dimensions({"k", "m", "n"}, {2, 2, 2});

  // Attempt to increase size to 256 units.
  ASSERT_DEATH(tensor = MKLTensor({"f", "g"}, {16, 16}), "");
}

// Checks that various invalid method arguments generate failures.
TEST(MKLTensorDeathTest, InvalidInput) {
  // Mismatched indices and dimensions.
  ASSERT_DEATH(MKLTensor({"a", "b", "c"}, {2, 2}), "");

  // Data vector size mismatch.
  std::vector<std::complex<float>> data(8);
  ASSERT_DEATH(MKLTensor({"a", "b"}, {2, 2}, data), "");

  MKLTensor tensor_abc({"a", "b", "c"}, {2, 2, 2});
  MKLTensor tensor_ac({"a", "c"}, {2, 2});

  // Projecting to index other than indices[0].
  ASSERT_DEATH(tensor_abc.project("b", 0, tensor_ac), "");

  // Projecting to bad index_value.
  ASSERT_DEATH(tensor_abc.project("a", 2, tensor_ac), "");

  // Projecting to too-small tensor.
  MKLTensor tensor_ac_small({"a", "c"}, {2, 1});
  ASSERT_DEATH(tensor_abc.project("a", 0, tensor_ac_small), "");

  // Renaming a non-existent index.
  ASSERT_DEATH(tensor_abc.rename_index("x", "y"), "");

  // Renaming an existing index to another existing index.
  ASSERT_DEATH(tensor_abc.rename_index("a", "b"), "");

  // Bundling on a partially-invalid set of indices.
  ASSERT_DEATH(tensor_abc.bundle({"a", "x"}, "ax"), "");

  // Bundling a valid but reordered set of indices.
  ASSERT_DEATH(tensor_abc.bundle({"b", "a"}, "ba"), "");

  // Reordering to too few indices.
  std::array<std::complex<float>, 256> scratch;
  ASSERT_DEATH(tensor_abc.reorder({"b", "a"}, scratch.data()), "");

  // Reordering to non-existent indices.
  ASSERT_DEATH(tensor_abc.reorder({"b", "y", "x"}, scratch.data()), "");

  MKLTensor tensor_cd({"c", "d"}, {2, 2});
  MKLTensor tensor_abd({"a", "b", "d"}, {2, 2, 2});

  // Reusing either tensor in multiplication.
  ASSERT_DEATH(multiply(tensor_abc, tensor_cd, tensor_abc, scratch.data()), "");
  ASSERT_DEATH(multiply(tensor_abc, tensor_cd, tensor_cd, scratch.data()), "");

  MKLTensor tensor_cd_large({"c", "d"}, {4, 4});

  // Mismatched index dimension in multiplication.
  ASSERT_DEATH(
      multiply(tensor_abc, tensor_cd_large, tensor_abd, scratch.data()), "");

  MKLTensor tensor_x({"x"}, {2});

  // Output tensor for multiplication is too small.
  ASSERT_DEATH(multiply(tensor_abc, tensor_cd, tensor_x, scratch.data()), "");
}

}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  return RUN_ALL_TESTS();
}
