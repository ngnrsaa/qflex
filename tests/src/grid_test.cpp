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
    EXPECT_THAT(msg, testing::HasSubstr("Cannot open grid file: invalid.txt."));
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
    EXPECT_THAT(msg, testing::HasSubstr("Grid size at line 2 is inconsistent "
                                        "with a width of 4 instead of 3."));
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
  std::vector<std::vector<std::size_t>> expected_off = {
      {0, 0}, {0, 3}, {2, 0}, {2, 2}, {2, 3}};
  EXPECT_EQ(grid.qubits_off, expected_off);
  EXPECT_EQ(grid.I, 3);
  EXPECT_EQ(grid.J, 4);
}

TEST(GridTest, ReadValidGrid6x2) {
  std::stringstream stream(kTestGrid_6x2);
  QflexGrid grid;
  grid.load(stream);
  std::vector<std::vector<std::size_t>> expected_off = {
      {0, 0}, {1, 1}, {4, 0}, {5, 0}, {5, 1}};
  EXPECT_EQ(grid.qubits_off, expected_off);
  EXPECT_EQ(grid.I, 6);
  EXPECT_EQ(grid.J, 2);
}

TEST(GridTest, ClearGridTest) {
  QflexGrid grid;
  grid.I = 4;
  grid.J = 8;
  grid.qubits_off = {{0, 0}, {1, 0}, {0, 1}};
  grid.clear();
  EXPECT_EQ(grid.I, 0);
  EXPECT_EQ(grid.J, 0);
  EXPECT_EQ(grid.qubits_off.size(), 0);
}

}  // namespace
}  // namespace qflex

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
