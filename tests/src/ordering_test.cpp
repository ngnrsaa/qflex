#include "ordering.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace qflex {
namespace {

TEST(OrderingExceptionTest, InvalidFilenameInput) {
  QflexOrdering ordering;
  std::string invalid_filename = "invalid.txt";
  try {
    ordering.load(invalid_filename);
  } catch (std::string msg) {
    EXPECT_THAT(msg,
                testing::HasSubstr("Cannot open ordering file: invalid.txt"));
  }
}

constexpr char kSimpleOrdering1[] = R"(# comment
cut () 0 1
# special character
expand *A 2
# tab instead of spaces
expand  A    4
# multiple spaces
expand      B   3
# trailing spaces
merge A B   
# spaces between parentheses
cut ( 0 ) 4
)";

constexpr char kSimpleOrdering2[] = R"(# comment
cut () 1 2
expand A 3
expand A 4
expand B 5
expand B 6
merge A B
)";

TEST(OrderingTest, LoadTest) {
  QflexOrdering ordering;

  // Check poorly formatted ordering.
  ordering.load(std::stringstream(kSimpleOrdering1));
  std::vector<std::string> check_ordering_1 = {"cut () 0 1", "expand A 2",
                                               "expand A 4", "expand B 3",
                                               "merge A B",  "cut (0) 4"};

  for (std::size_t i = 0; i < std::size(check_ordering_1); ++i) {
    ASSERT_EQ(check_ordering_1[i], ordering.instructions[i]);
  }

  // Check normal ordering and that load() clears previous ordering.
  ordering.load(std::stringstream(kSimpleOrdering2));
  std::vector<std::string> check_ordering_2 = {"cut () 1 2", "expand A 3",
                                               "expand A 4", "expand B 5",
                                               "expand B 6", "merge A B"};
  for (std::size_t i = 0; i < std::size(check_ordering_2); ++i) {
    ASSERT_EQ(check_ordering_2[i], ordering.instructions[i]);
  }
}

}  // namespace
}  // namespace qflex

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
