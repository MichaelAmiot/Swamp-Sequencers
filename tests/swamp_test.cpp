#include "../include/SwampSeqLib/suffix_array.h"
#include "SwampSeqLib/genome_mapper.h"
#include "SwampSeqLib/suffix_tree.h"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <random>

void sortResults(std::vector<STSearchResult> &res) {
  std::sort(res.begin(), res.end(),
            [](const STSearchResult &a, const STSearchResult &b) {
              return a.offset < b.offset;
            });
}

const std::vector<char> possibleChars{'A', 'C', 'T', 'G', 'N'};

TEST(AlwaysTrue, AlwaysTrue) { ASSERT_EQ(0, 0); }

TEST(SuffixArray, BasicString) {
  const std::string text = "banana";
  SuffixArray SA(text);

  std::vector<size_t> expected{5, 3, 1, 0, 4, 2};
  EXPECT_EQ(SA.sa(), expected);
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

TEST(SuffixTree, BasicString) {
  const std::string text = "banana";
  SuffixTree ST(text);

  // Searching for "a" to verify multiple leaf discovery
  std::vector<STSearchResult> res = ST.search("a");
  sortResults(res);

  ASSERT_EQ(res.size(), 3);
  EXPECT_EQ(res[0].offset, 1);
  EXPECT_EQ(res[1].offset, 3);
  EXPECT_EQ(res[2].offset, 5);
}

TEST(SuffixTree, GenomeString) {
  const std::string genome = "ACGGTCGATCACGAT";
  SuffixTree ST(genome);

  std::vector<STSearchResult> expected{{0, 3}, {10, 3}};
  auto res = ST.search("ACG");
  sortResults(res);

  ASSERT_EQ(res.size(), expected.size());
  for (size_t i = 0; i < res.size(); ++i) {
    EXPECT_EQ(res[i].offset, expected[i].offset);
    EXPECT_EQ(res[i].length, expected[i].length);
  }
}

TEST(SuffixTree, Repetition) {
  std::string genome(1000, 'A');
  SuffixTree ST(genome);

  auto res = ST.search("AA");
  sortResults(res);

  // "AA" appears at every index from 0 to 998
  ASSERT_EQ(res.size(), 999);
  for (int64_t i = 0; i < 999; i++) {
    EXPECT_EQ(res[i].offset, i);
    EXPECT_EQ(res[i].length, 2);
  }
}

TEST(SuffixTree, Overlapping) {
  std::string text = "ba";
  for (size_t i = 0; i < 100; i++)
    text += "na";

  SuffixTree ST(text);

  // Pattern "ana" occurs starting at indices 1, 3, 5...
  auto res = ST.search("ana");
  sortResults(res);

  size_t expected_count = (text.size() - 1) / 2; // "ana" in "bananana..."

  for (size_t i = 0; i < res.size(); i++) {
    EXPECT_EQ(res[i].offset, (int64_t)(i * 2 + 1));
    EXPECT_EQ(res[i].length, 3);
  }
}

TEST(SuffixTree, NoMatch) {
  std::string genome(1000, 'A');
  SuffixTree ST(genome);

  auto res = ST.search("NO MATCH");
  EXPECT_EQ(res.size(), 0);
}

TEST(SuffixTree, Boundaries) {
  const std::string text = "GATTACA";
  SuffixTree ST(text);

  auto resStart = ST.search("GAT");
  ASSERT_EQ(resStart.size(), 1);
  EXPECT_EQ(resStart[0].offset, 0);

  auto resEnd = ST.search("ACA");
  ASSERT_EQ(resEnd.size(), 1);
  EXPECT_EQ(resEnd[0].offset, 4);
}

TEST(SuffixTree, AllMatch) {
  std::string genome(100, 'A');
  SuffixTree ST(genome);

  auto res = ST.search("A");
  EXPECT_EQ(res.size(), 100);
}

TEST(SuffixTree, EmptyInput) {
  EXPECT_NO_THROW({
    SuffixTree ST("");
    auto res = ST.search("A");
    EXPECT_EQ(res.size(), 0);
    EXPECT_EQ(ST.size(), 0);
  });
}

TEST(SuffixTree, SingleChar) {
  SuffixTree ST("G");
  auto resG = ST.search("G");
  EXPECT_EQ(resG.size(), 1);
  if (!resG.empty())
    EXPECT_EQ(resG[0].offset, 0);

  EXPECT_EQ(ST.search("C").size(), 0);
}

TEST(SuffixTree, InternalState) {
  SuffixTree ST("banana");
  EXPECT_TRUE(ST.ready());
  EXPECT_EQ(ST.size(), 6);
}
