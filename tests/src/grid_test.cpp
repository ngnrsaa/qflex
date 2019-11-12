#include "grid.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace qflex {
namespace {

TEST(GridLoadTest, InvalidFilenameInput) {
    QflexGrid grid; 
    std::string invalid_filename = "invalid.txt";

    EXPECT_ANY_THROW(grid.load(invalid_filename));

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

TEST(ReadGridTest, ValidGrid3x4) {
  std::stringstream stream(kTestGrid_3x4);
  QflexGrid grid;
  grid.load(stream);
  std::vector<std::vector<int>> expected_off = {
      {0, 0}, {0, 3}, {2, 0}, {2, 2}, {2, 3}};
  EXPECT_EQ(grid.qubits_off, expected_off);
}

TEST(ReadGridTest, ValidGrid6x2) {
  std::stringstream stream(kTestGrid_6x2);
  QflexGrid grid;
  grid.load(stream);
  std::vector<std::vector<int>> expected_off = {
      {0, 0}, {1, 1}, {4, 0}, {5, 0}, {5, 1}};
  EXPECT_EQ(grid.qubits_off, expected_off);
}


}  // namespace
}  // namespace qflex

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}