#include "tensor.h"

#include <array>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace qflex {
namespace {

using ::testing::Eq;
using ::testing::Pointwise;

/**
 * Simple reordering routine to swap axes of a tensor
 * @param array input array
 * @param initial_indices set of tags corresponding to the indices of the input
 * tensor
 * @param final_indices set of tags corresponding to the desired re-ordering
 * @return a new tensor with the desired order of indices
 */
template <typename array_type, typename initial_indices_type,
          typename final_indices_type>
array_type simple_reordering(const array_type &array,
                             const initial_indices_type &initial_indices,
                             const final_indices_type &final_indices) {
  // Find map
  std::vector<std::size_t> map;
  for (std::size_t i = 0; i < std::size(initial_indices); ++i)
    map.push_back(
        std::size(initial_indices) -
        std::distance(std::begin(initial_indices),
                      std::find(std::begin(initial_indices),
                                std::end(initial_indices), final_indices[i])) -
        1);
  std::reverse(std::begin(map), std::end(map));

  array_type out(std::size(array));
  for (std::size_t i = 0; i < std::size(array); ++i) {
    // Apply permutation to indices
    std::size_t j = 0;
    for (std::size_t p1 = 0; p1 < std::size(map); ++p1) {
      std::size_t p2 = map[p1];
      j ^= ((i >> p1) & std::size_t(1)) << p2;
    }

    // Assign
    out[i] = array[j];
  }

  return out;
}

// Creates an empty tensor and runs basic sanity checks on it.
TEST(TensorTest, EmptyTensor) {
  std::vector<std::string> indices = {"a", "b"};
  std::vector<size_t> dimensions = {2, 4};

  // Automatically allocates new space.
  Tensor tensor(indices, dimensions);
  ASSERT_EQ(tensor.get_indices(), indices);
  ASSERT_EQ(tensor.get_dimensions(), dimensions);
  ASSERT_EQ(tensor.size(), 8ul);
  ASSERT_EQ(tensor.num_zeros(), 8ul);

  // Verify that data is initialized to zero.
  std::vector<std::complex<float>> read_data(tensor.data(),
                                             tensor.data() + tensor.size());
  for (std::size_t i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), 0);
    ASSERT_FLOAT_EQ(read_data[i].imag(), 0);
  }

  const auto dict = tensor.get_index_to_dimension();
  ASSERT_EQ(dict.at("a"), 2ul);
  ASSERT_EQ(dict.at("b"), 4ul);

  // Test index renaming as well.
  tensor.rename_index("a", "c");
  std::vector<std::string> expected_indices = {"c", "b"};
  ASSERT_EQ(tensor.get_indices(), expected_indices);
}

// Loads a tensor from data and runs basic sanity checks on it.
TEST(TensorTest, LoadData) {
  std::vector<std::string> indices = {"a", "b"};
  std::vector<size_t> dimensions = {2, 2};
  std::vector<std::complex<float>> data = {
      std::complex<float>(0, 1),
      std::complex<float>(2, 3),
      std::complex<float>(4, 5),
      std::complex<float>(6, 7),
  };
  Tensor tensor(indices, dimensions, data);
  ASSERT_EQ(tensor.get_indices(), indices);
  ASSERT_EQ(tensor.get_dimensions(), dimensions);
  ASSERT_EQ(tensor.size(), 4ul);
  ASSERT_EQ(tensor.num_zeros(), 0ul);

  // Verify that data is read in from vector as expected.
  std::vector<std::complex<float>> read_data(tensor.data(),
                                             tensor.data() + tensor.size());
  for (std::size_t i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), data[i].real());
    ASSERT_FLOAT_EQ(read_data[i].imag(), data[i].imag());
  }

  // Test scalar multiplication as well.
  tensor.scalar_multiply(std::complex<float>(0, 1));
  read_data = std::vector<std::complex<float>>(tensor.data(),
                                               tensor.data() + tensor.size());
  for (std::size_t i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), -1 * data[i].imag());
    ASSERT_FLOAT_EQ(read_data[i].imag(), data[i].real());
  }
}

// Projects a tensor onto a single value of an index and verifies that the
// output tensor only contains data from that slice of the original tensor.
TEST(TensorTest, TensorProjection) {
  std::vector<std::string> indices = {"a", "b", "c"};
  std::vector<size_t> dimensions = {2, 2, 2};
  std::vector<std::complex<float>> data;
  for (std::size_t i = 0; i < 8; i++) {
    data.push_back(std::complex<float>(i, 0));
  }

  Tensor tensor(indices, dimensions, data);
  std::vector<std::string> expected_indices = {"b", "c"};
  std::vector<size_t> expected_dimensions = {2, 2};

  Tensor projection_tensor_1({"x", "y"}, {2, 2});
  tensor.project("a", 1, projection_tensor_1);
  ASSERT_EQ(projection_tensor_1.get_indices(), expected_indices);
  ASSERT_EQ(projection_tensor_1.get_dimensions(), expected_dimensions);

  Tensor projection_tensor_2({""}, {64});
  tensor.project("a", 1, projection_tensor_2);
  ASSERT_EQ(projection_tensor_2.get_indices(), expected_indices);
  ASSERT_EQ(projection_tensor_2.get_dimensions(), expected_dimensions);

  const std::size_t psize = projection_tensor_1.size();

  const std::vector<std::complex<float>> proj_data_1(
      projection_tensor_1.data(), projection_tensor_1.data() + psize);
  const std::vector<std::complex<float>> proj_data_2(
      projection_tensor_2.data(), projection_tensor_2.data() + psize);
  for (std::size_t i = 0; i < psize; ++i) {
    ASSERT_FLOAT_EQ(proj_data_1[i].real(), data[i + psize].real());
    ASSERT_FLOAT_EQ(proj_data_2[i].real(), data[i + psize].real());
  }
}

// Bundles indices of a tensor and verifies that dimensions change while data
// remains unaffected.
TEST(TensorTest, IndexBundling) {
  std::vector<std::string> indices = {"a", "b", "c", "d"};
  std::vector<size_t> dimensions = {2, 2, 2, 2};
  std::vector<std::complex<float>> data;
  for (std::size_t i = 0; i < 16; i++) {
    data.push_back(std::complex<float>(i, 0));
  }

  Tensor tensor(indices, dimensions, data);
  tensor.bundle({"a", "b", "c"}, "abc");
  std::vector<std::string> expected_indices = {"abc", "d"};
  std::vector<size_t> expected_dimensions = {8, 2};
  ASSERT_EQ(tensor.get_indices(), expected_indices);
  ASSERT_EQ(tensor.get_dimensions(), expected_dimensions);

  std::vector<std::complex<float>> read_data(tensor.data(),
                                             tensor.data() + tensor.size());
  for (std::size_t i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), data[i].real());
  }
}

// Reorders indices of a tensor and verifies that data changes accordingly.
TEST(TensorTest, SimpleIndexReordering) {
  std::vector<std::string> indices = {"a", "b", "c"};
  std::vector<size_t> dimensions = {2, 2, 2};
  std::vector<std::complex<float>> data;
  for (std::size_t i = 0; i < 8; i++) {
    data.push_back(std::complex<float>(i, 0));
  }

  Tensor tensor(indices, dimensions, data);
  std::vector<std::string> expected_indices = {"b", "c", "a"};
  std::array<std::complex<float>, 8> scratch;
  try {
    tensor.reorder(expected_indices, scratch.data());
  } catch (std::string msg) {
    FAIL()
        << "Expected tensor.reorder() to succeed but failed with error msg:  "
        << msg << std::endl;
  }
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
  for (std::size_t i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), expected_data[i].real());
  }
}

// Tests a reordering of ten indices needing a single right move.
// <empty> | abcde | fghij
// -> <empty> | jegcf | bhida
TEST(TensorTest, RightTenIndicesReordering) {
  std::vector<std::string> indices = {"a", "b", "c", "d", "e",
                                      "f", "g", "h", "i", "j"};
  std::vector<size_t> dimensions = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  std::vector<std::complex<float>> data;
  for (int i = 0; i < 1024; i++) {
    data.push_back(std::complex<float>(i, 0));
  }

  Tensor tensor(indices, dimensions, data);
  std::vector<std::string> expected_indices = {"j", "e", "g", "c", "f",
                                               "b", "h", "i", "d", "a"};
  std::array<std::complex<float>, 1024> scratch;
  try {
    tensor.reorder(expected_indices, scratch.data());
  } catch (std::string msg) {
    FAIL()
        << "Expected tensor.reorder() to succeed but failed with error msg:  "
        << msg << std::endl;
  }
  ASSERT_EQ(tensor.get_indices(), expected_indices);
  ASSERT_EQ(tensor.get_dimensions(), dimensions);

  // Check Tensor data.
  data = simple_reordering(data, indices, expected_indices);

  for (std::size_t i = 0; i < std::size(data); ++i) {
    ASSERT_EQ(data[i], tensor.data()[i]);
  }
}

// Tests a reordering needing a single right move.
// ab | cdefg | hijkl
// -> ab | ldefh | gijkc
TEST(TensorTest, RightTwelveIndicesReordering) {
  std::vector<std::string> indices = {"a", "b", "c", "d", "e", "f",
                                      "g", "h", "i", "j", "k", "l"};
  std::vector<size_t> dimensions = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  std::vector<std::complex<float>> data;
  for (int i = 0; i < 4096; i++) {
    data.push_back(std::complex<float>(i, 0));
  }
  Tensor tensor(indices, dimensions, data);
  std::vector<std::string> expected_indices = {"a", "b", "l", "d", "e", "f",
                                               "h", "g", "i", "j", "k", "c"};
  std::array<std::complex<float>, 4096> scratch;
  try {
    tensor.reorder(expected_indices, scratch.data());
  } catch (std::string msg) {
    FAIL()
        << "Expected tensor.reorder() to succeed but failed with error msg:  "
        << msg << std::endl;
  }
  ASSERT_EQ(tensor.get_indices(), expected_indices);
  ASSERT_EQ(tensor.get_dimensions(), dimensions);

  // Check Tensor data.
  data = simple_reordering(data, indices, expected_indices);

  for (std::size_t i = 0; i < std::size(data); ++i) {
    ASSERT_EQ(data[i], tensor.data()[i]);
  }
}

// Tests a reordering needing a single left move.
// ab | cdefg | hijkl
// -> ga | cdefb | hijkl
TEST(TensorTest, LeftTwelveIndicesReordering) {
  std::vector<std::string> indices = {"a", "b", "c", "d", "e", "f",
                                      "g", "h", "i", "j", "k", "l"};
  std::vector<size_t> dimensions = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  std::vector<std::complex<float>> data;
  for (int i = 0; i < 4096; i++) {
    data.push_back(std::complex<float>(i, 0));
  }
  Tensor tensor(indices, dimensions, data);
  std::vector<std::string> expected_indices = {"g", "a", "c", "d", "e", "f",
                                               "b", "h", "i", "j", "k", "l"};
  std::array<std::complex<float>, 4096> scratch;
  try {
    tensor.reorder(expected_indices, scratch.data());
  } catch (std::string msg) {
    FAIL()
        << "Expected tensor.reorder() to succeed but failed with error msg:  "
        << msg << std::endl;
  }
  ASSERT_EQ(tensor.get_indices(), expected_indices);
  ASSERT_EQ(tensor.get_dimensions(), dimensions);

  // Check Tensor data.
  data = simple_reordering(data, indices, expected_indices);

  for (std::size_t i = 0; i < std::size(data); ++i) {
    ASSERT_EQ(data[i], tensor.data()[i]);
  }
}

// clang-format off
// Tests a reordering needing a single left move. Even though index 'h' should be
// in the rightmost grouping, left_reorder() will check and see if it can perform
// the reordering by making a left move with up to the 8th index from the left. 
// abcdefgh | ijkl
// -> hbcdefga | ijkl
// clang-format on
TEST(TensorTest, SlowLeftTwelveIndicesReordering) {
  std::vector<std::string> indices = {"a", "b", "c", "d", "e", "f",
                                      "g", "h", "i", "j", "k", "l"};
  std::vector<size_t> dimensions = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  std::vector<std::complex<float>> data;
  for (int i = 0; i < 4096; i++) {
    data.push_back(std::complex<float>(i, 0));
  }
  Tensor tensor(indices, dimensions, data);
  std::vector<std::string> expected_indices = {"h", "b", "c", "d", "e", "f",
                                               "g", "a", "i", "j", "k", "l"};
  std::array<std::complex<float>, 4096> scratch;
  try {
    tensor.reorder(expected_indices, scratch.data());
  } catch (std::string msg) {
    FAIL()
        << "Expected tensor.reorder() to succeed but failed with error msg:  "
        << msg << std::endl;
  }
  ASSERT_EQ(tensor.get_indices(), expected_indices);
  ASSERT_EQ(tensor.get_dimensions(), dimensions);

  // Check Tensor data.
  data = simple_reordering(data, indices, expected_indices);

  for (std::size_t i = 0; i < std::size(data); ++i) {
    ASSERT_EQ(data[i], tensor.data()[i]);
  }
}

// Tests a reordering needing a left and a right move.
// Left: ab | cdefg | hijkl -> cd | abefg | hijkl
// Right: cd | abefg | hijkl -> cd | hijkl | efgab
TEST(TensorTest, LeftRightIndexReordering) {
  std::vector<std::string> indices = {"a", "b", "c", "d", "e", "f",
                                      "g", "h", "i", "j", "k", "l"};
  std::vector<size_t> dimensions = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  std::vector<std::complex<float>> data;
  for (int i = 0; i < 4096; i++) {
    data.push_back(std::complex<float>(i, 0));
  }
  Tensor tensor(indices, dimensions, data);
  std::vector<std::string> expected_indices = {"c", "d", "h", "i", "j", "k",
                                               "l", "e", "f", "g", "a", "b"};
  std::array<std::complex<float>, 4096> scratch;
  try {
    tensor.reorder(expected_indices, scratch.data());
  } catch (std::string msg) {
    FAIL()
        << "Expected tensor.reorder() to succeed but failed with error msg:  "
        << msg << std::endl;
  }

  // Check Tensor data.
  data = simple_reordering(data, indices, expected_indices);

  for (std::size_t i = 0; i < std::size(data); ++i) {
    ASSERT_EQ(data[i], tensor.data()[i]);
  }
}

// Tests a worse case index reordering needing a left, right, and left move.
// Left: ab | cdefg | hijkl -> ac | dfgbe | hijkl
// Right: ac | dfgbe | hijkl -> ac | dfgkl | hbeij
// Left: ac | dfgkl | hbeij -> kc | aldgf | hbeij
TEST(TensorTest, WorstCaseIndexReordering) {
  std::vector<std::string> indices = {"a", "b", "c", "d", "e", "f",
                                      "g", "h", "i", "j", "k", "l"};
  std::vector<size_t> dimensions = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  std::vector<std::complex<float>> data;
  for (int i = 0; i < 4096; i++) {
    data.push_back(std::complex<float>(i, 0));
  }

  Tensor tensor(indices, dimensions, data);
  std::vector<std::string> expected_indices = {"k", "c", "a", "l", "d", "g",
                                               "f", "h", "b", "e", "i", "j"};
  std::array<std::complex<float>, 4096> scratch;
  try {
    tensor.reorder(expected_indices, scratch.data());
  } catch (std::string msg) {
    FAIL()
        << "Expected tensor.reorder() to succeed but failed with error msg:  "
        << msg << std::endl;
  }
  ASSERT_EQ(tensor.get_indices(), expected_indices);
  ASSERT_EQ(tensor.get_dimensions(), dimensions);

  // Check Tensor data.
  data = simple_reordering(data, indices, expected_indices);

  for (std::size_t i = 0; i < std::size(data); ++i) {
    ASSERT_EQ(data[i], tensor.data()[i]);
  }
}

// Multiplies two tensors and verify shape, indices, and data of the result.
TEST(TensorTest, Multiply) {
  std::vector<std::string> indices_a = {"a", "b", "c"};
  std::vector<size_t> dimensions_a = {2, 2, 2};
  std::vector<std::complex<float>> data_a;
  for (std::size_t i = 0; i < 8; i++) {
    data_a.push_back(std::complex<float>(i, 0));
  }
  Tensor tensor_a(indices_a, dimensions_a, data_a);

  std::vector<std::string> indices_b = {"b", "c", "d"};
  std::vector<size_t> dimensions_b = {2, 2, 2};
  std::vector<std::complex<float>> data_b;
  for (std::size_t i = 8; i > 0; i--) {
    data_b.push_back(std::complex<float>(i, 0));
  }
  Tensor tensor_b(indices_b, dimensions_b, data_b);

  std::vector<std::string> indices_c = {"x"};
  std::vector<size_t> dimensions_c = {16};
  Tensor tensor_c(indices_c, dimensions_c);

  std::array<std::complex<float>, 16> scratch;
  multiply(tensor_a, tensor_b, tensor_c, scratch.data());
  std::vector<std::string> expected_indices = {"a", "d"};
  std::vector<size_t> expected_dimensions = {2, 2};
  ASSERT_EQ(tensor_c.get_indices(), expected_indices);
  ASSERT_EQ(tensor_c.get_dimensions(), expected_dimensions);

  std::vector<std::complex<float>> expected_data = {20, 14, 100, 78};
  std::vector<std::complex<float>> read_data(tensor_c.data(),
                                             tensor_c.data() + tensor_c.size());
  for (std::size_t i = 0; i < read_data.size(); ++i) {
    ASSERT_FLOAT_EQ(read_data[i].real(), expected_data[i].real());
  }
}

// Verifies that a tensor retains its initialized capacity.
TEST(TensorExceptionTest, Capacity) {
  std::vector<std::string> indices = {"a", "b"};
  std::vector<size_t> dimensions = {2, 4};

  // Create an uninitialized tensor.
  Tensor tensor;

  // Allocate 64 units of space, which should always be available.
  tensor = Tensor({""}, {64});

  // Reduce size to 16 units; capacity is still 64 units.
  tensor.set_indices_and_dimensions({"a", "b"}, {4, 4});

  // Increase size to 32 units; capacity is still 64 units.
  tensor.set_indices_and_dimensions({"x"}, {32});

  // Reduce size to 2 units; capacity is still 64 units.
  tensor.set_indices_and_dimensions({"z"}, {2});

  // Increase size to 8 units; capacity is still 64 units.
  tensor.set_indices_and_dimensions({"k", "m", "n"}, {2, 2, 2});

  // Moving tensors is always allowed
  tensor = Tensor({"f", "g"}, {16, 16});

  // Attempt to increase size to 4096 units.
  try {
    const auto new_tensor = Tensor({"f", "g", "h"}, {16, 16, 16});
    tensor = new_tensor;
    FAIL() << "Expected Tensor() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr(
                         "The total allocated space: 256, is insufficient for "
                         "the requested tensor dimensions: 4096."));
  }
}

// Checks that various invalid method arguments generate failures.
TEST(TensorExceptionTest, InvalidInput) {
  // Mismatched indices and dimensions.
  try {
    Tensor({"a", "b", "c"}, {2, 2});
    FAIL() << "Expected Tensor() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("The number of indices: 3, and number "
                                        "of dimensions: 2, should be equal."));
  }

  // Data vector size mismatch.
  std::vector<std::complex<float>> data(8);
  try {
    Tensor({"a", "b"}, {2, 2}, data);
    FAIL() << "Expected Tensor() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg,
        testing::HasSubstr(
            "The vector data size: 8, has to match the size of the Tensor: 4"));
  }

  Tensor tensor_abc({"a", "b", "c"}, {2, 2, 2});
  Tensor tensor_ac({"a", "c"}, {2, 2});

  // Projecting to index other than indices[0].
  try {
    tensor_abc.project("b", 0, tensor_ac);
    FAIL() << "Expected project() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr(
                         "Index: 'b' has to be equal to indices[0]: 'a'."));
  }

  // Projecting to bad index_value.
  try {
    tensor_abc.project("a", 2, tensor_ac);
    FAIL() << "Expected project() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg, testing::HasSubstr("index_value: 2 must be contained in [0, 2)."));
  }

  // Projecting to too-small tensor.
  Tensor tensor_ac_small({"a", "c"}, {2, 1});
  try {
    tensor_abc.project("a", 0, tensor_ac_small);
    FAIL() << "Expected project() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg, testing::HasSubstr("The total allocated space: 2, is insufficient "
                                "for the requested tensor dimensions: 4."));
  }

  // Renaming a non-existent index.
  try {
    tensor_abc.rename_index("x", "y");
    FAIL() << "Expected rename_index() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg,
                testing::HasSubstr("old_name: x, has to be a valid index."));
  }

  // Renaming an existing index to another existing index.
  try {
    tensor_abc.rename_index("a", "b");
    FAIL() << "Expected rename_index() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg, testing::HasSubstr("new_name: b, cannot be an existing index."));
  }

  // Bundling on a partially-invalid set of indices.
  try {
    tensor_abc.bundle({"a", "x"}, "ax");
    FAIL() << "Expected bundle() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("indices_to_bundle: {a, x} has to be "
                                        "contained in indices: {a, b, c}."));
  }

  // Bundling a valid but reordered set of indices.
  try {
    tensor_abc.bundle({"b", "a"}, "ba");
    FAIL() << "Expected bundle() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("indices_to_bundle: {b, a} must be in "
                                        "its original order: {a, b}."));
  }

  // Reordering to too few indices.
  std::array<std::complex<float>, 256> scratch;
  try {
    tensor_abc.reorder({"b", "a"}, scratch.data());
    FAIL() << "Expected reorder() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg,
                testing::HasSubstr("new_ordering: {b, a} must be a reordering "
                                   "of current indices: {a, b, c}."));
  }

  // Reordering to non-existent indices.
  try {
    tensor_abc.reorder({"b", "y", "x"}, scratch.data());
    FAIL() << "Expected reorder() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg, testing::HasSubstr("new_ordering: {b, y, x} must be a reordering "
                                "of current indices: {a, b, c}."));
  }

  // Scratch copy passed to reordering cannot be null pointer.
  try {
    tensor_abc.reorder({"a", "b"}, nullptr);
    FAIL() << "Expected reorder() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("Scratch copy must be non-null."));
  }

  Tensor tensor_cd({"c", "d"}, {2, 2});
  Tensor tensor_abd({"a", "b", "d"}, {2, 2, 2});

  // Reusing either tensor in multiplication.
  try {
    multiply(tensor_abc, tensor_cd, tensor_abc, scratch.data());
    FAIL() << "Expected multiply() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg,
                testing::HasSubstr("A and C cannot be the same tensor: Tensor "
                                   "of rank 3: a -> 2, b -> 2, c -> 2"));
  }
  try {
    multiply(tensor_abc, tensor_cd, tensor_cd, scratch.data());
    FAIL() << "Expected multiply() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("B and C cannot be the same tensor: "
                                        "Tensor of rank 2: c -> 2, d -> 2"));
  }

  Tensor tensor_cd_large({"c", "d"}, {4, 4});

  // Mismatched index dimension in multiplication.
  try {
    multiply(tensor_abc, tensor_cd_large, tensor_abd, scratch.data());
    FAIL() << "Expected multiply() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr(
                         "Common indices must have matching dimensions, but at "
                         "index: 0, the dimensions are A: 2, B: 4."));
  }

  Tensor tensor_x({"x"}, {2});

  // Output tensor for multiplication is too small.
  try {
    multiply(tensor_abc, tensor_cd, tensor_x, scratch.data());
    FAIL() << "Expected multiply() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(
        msg, testing::HasSubstr(
                 "C: 2 doesn't have enough space for the product of A*B: 8."));
  }

  // Scratch copy passed to multiply cannot be null pointer.
  try {
    multiply(tensor_abc, tensor_cd, tensor_x, nullptr);
    FAIL() << "Expected multiply() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("Scratch copy must be non-null."));
  }
}

// Testing this function by direct call because too nested to test by calling
// Tensor::reorder
TEST(TensorExceptionTest, GenerateBinaryReorderingMapInvalidInput) {
  const std::vector<std::size_t> map_old_to_new_idxpos = {1, 2};
  std::vector<std::size_t> map_old_to_new_position = {1, 2, 3};

  // Size of map must be equal to 2 ^ (number of indices).
  try {
    _generate_binary_reordering_map(map_old_to_new_idxpos,
                                    map_old_to_new_position);
    FAIL()
        << "Expected _generate_binary_reordering_map() to throw an exception.";
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr(
                         "Size of map: 3 must be equal to 2^num_indices: 4."));
  }
}

}  // namespace
}  // namespace qflex

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  return RUN_ALL_TESTS();
}
