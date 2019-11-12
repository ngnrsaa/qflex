#include "grid.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace qflex {
namespace {

TEST(GridExceptionTest, InvalidFilenameInput) {
  QflexGrid grid;
  std::string invalid_filename = "invalid.txt";
  try {
    grid.load(invalid_filename);
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("Cannot open grid file: invalid.txt"));
  }
}

// Invalid grid layout.
constexpr char kInvalidTestGrid[] = R"(1 1 1
                                       1 1 1 1
                                       1 1 1)";

TEST(GridExceptionTest, InvalidGridFormat) {
  std::stringstream stream(kInvalidTestGrid);
  QflexGrid grid;
  try {
    grid.load(stream);
  } catch (std::string msg) {
    EXPECT_THAT(msg, testing::HasSubstr("Grid size is inconsistent."));
  }
}

// Grid layout with trailing whitespace.
constexpr char kTestGrid_3x4[] = R"(0 1 1 0
                                    1 1 1 1
                                    0 1 0 0)";

constexpr char kTestGrid_6x2[] = R"(0 1
                                    1 0
                                    1 1
                                    1 1
                                    0 1
                                    0 0)";

TEST(GridTest, ReadValidGrid3x4) {
  std::stringstream stream(kTestGrid_3x4);
  QflexGrid grid;
  grid.load(stream);
  std::vector<std::vector<int>> expected_off = {
      {0, 0}, {0, 3}, {2, 0}, {2, 2}, {2, 3}};
  EXPECT_EQ(grid.qubits_off, expected_off);
  EXPECT_EQ(grid.I, 3);
  EXPECT_EQ(grid.J, 4);
}

TEST(GridTest, ReadValidGrid6x2) {
  std::stringstream stream(kTestGrid_6x2);
  QflexGrid grid;
  grid.load(stream);
  std::vector<std::vector<int>> expected_off = {
      {0, 0}, {1, 1}, {4, 0}, {5, 0}, {5, 1}};
  EXPECT_EQ(grid.qubits_off, expected_off);
  EXPECT_EQ(grid.I, 6);
  EXPECT_EQ(grid.J, 2);
}

}  // namespace
}  // namespace qflex

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}