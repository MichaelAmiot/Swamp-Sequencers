#include "SwampSeqLib/genome_mapper.h"
#include "../include/SwampSeqLib/suffix_array.h"
#include "SwampSeqLib/suffix_tree.h"
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <random>

const std::vector<char> possibleChars{'A', 'C', 'T', 'G', 'N'};

TEST(AlwaysTrue, AlwaysTrue) { ASSERT_EQ(0, 0); }

TEST(SuffixArray, BasicString) {
  const std::string text = "banana";
  SuffixArray SA(text);

  std::vector<size_t> expected{5, 3, 1, 0, 4, 2};
  EXPECT_EQ(SA.getArray(), expected);
}

TEST(SuffixArray, GenomeString) {
  const std::string genome = "ACGGTCGATCACGAT";
  SuffixArray SA(genome);

  std::vector<SearchResult> expected{{0, 3}, {10, 3}};
  EXPECT_EQ(SA.search("ACG"), expected);
}

TEST(SuffixArray, Repetition) {
  std::string genome = "";
  for (size_t i = 0; i < 1000; i++)
    genome += 'A';
  SuffixArray SA(genome);

  std::vector<SearchResult> expected;
  for (size_t i = 0; i < 999; i += 1)
    expected.push_back({i, 2});
  auto res = SA.search("AA");
  for (size_t i = 0; i < 500; i++)
    EXPECT_EQ(res[i].offset, expected[i].offset);
}

TEST(SuffixArray, Overlapping) {
  std::string text = "ba";
  for (size_t i = 0; i < 100; i++)
    text += "na";

  SuffixArray SA(text);

  std::vector<SearchResult> expected;
  for (size_t i = 0; i < (text.size() - 1) / 3; i++) {
    size_t offset = i * 2 + 1;
    expected.push_back({offset, 3});
  }

  auto res = SA.search("ana");
  for (size_t i = 0; i < (text.size() - 1) / 3; i++) {
    EXPECT_EQ(res[i].offset, expected[i].offset);
  }
}

TEST(SuffixArray, NoMatch) {
  std::string genome = "";
  for (size_t i = 0; i < 1000; i++)
    genome += "A";

  SuffixArray SA(genome);
  auto res = SA.search("NO MATCH");

  EXPECT_EQ(res.size(), 0);
}

TEST(SuffixArray, Boundaries) {
  const std::string text = "GATTACA";
  SuffixArray SA(text);

  auto resStart = SA.search("GAT");
  ASSERT_EQ(resStart.size(), 1);
  EXPECT_EQ(resStart[0].offset, 0);

  auto resEnd = SA.search("ACA");
  ASSERT_EQ(resEnd.size(), 1);
  EXPECT_EQ(resEnd[0].offset, 4);
}

TEST(SuffixArray, AllMatch) {
  std::string genome(100, 'A');
  SuffixArray SA(genome);

  auto res = SA.search("A");
  EXPECT_EQ(res.size(), 100);
}

// Constructor shouldn't segfault on empty input
TEST(SuffixArray, EmptyInput) {
  EXPECT_NO_THROW({
    SuffixArray SA("");
    auto res = SA.search("A");
    EXPECT_EQ(res.size(), 0);
  });
}

TEST(SuffixArray, SingleChar) {
  SuffixArray SA("G");
  EXPECT_EQ(SA.search("G").size(), 1);
  EXPECT_EQ(SA.search("C").size(), 0);
}
