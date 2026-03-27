#include "SwampSeqLib/suffix_array.h"
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <iostream>

// SearchResult comparison for testing
bool SearchResult::operator==(const SearchResult &other) const {
  return ((this->length == other.length) && (this->offset == other.offset));
}

// ─────────────────────────────────────────────────────────────────────────────
// Construction
// ─────────────────────────────────────────────────────────────────────────────

SuffixArray::SuffixArray(const GenomeMapper &mapper) {
  if (!mapper.isValid())
    throw std::runtime_error("SuffixArray: GenomeMapper is not valid.");

  _data = mapper.data();
  _num = mapper.size();

  buildSuffixArray();
  _ready = true;
}

SuffixArray::SuffixArray(const std::string &text) {
  _owned = text; // take a copy so lifetime is self-managed
  _data = _owned.c_str();
  _num = _owned.size();

  buildSuffixArray();
  _ready = true;
}

// ─────────────────────────────────────────────────────────────────────────────
// SA-IS  —  O(n) time, O(n) space
//
// Reference: Nong, Zhang & Chan (2009)
//   "Linear Suffix Array Construction by Almost Pure Induced Sorting"
//
// High-level outline
// ──────────────────
//   Classify every position as S-type or L-type.
//   Locate Left-Most-S (LMS) substrings.
//   Induced-sort the LMS suffixes into bucket heads/tails to get an
//     approximate SA.
//   Compact the sorted LMS substrings into a reduced string s1 and recurse.
//   Use the exact LMS order from the recursion to re-run induced sort on
//     the original string, yielding the final SA.
//
// Alphabet convention
// ───────────────────
//   The sentinel character 0 is implicitly appended (it never appears in the
//    real text and is lexicographically smallest).
//   The function works entirely on integer arrays so it handles both the
//    byte-alphabet first call and the reduced-alphabet recursive calls.
// ─────────────────────────────────────────────────────────────────────────────

// ── Helpers (file-scope) ─────────────────────────────────

// Return bucket sizes: buckets[c] = number of occurrences of symbol c.
std::vector<size_t> getBuckets(const std::vector<int64_t> &symbols,
                               int64_t alphabetSize) {
  const size_t bucketCount = static_cast<size_t>(alphabetSize);
  std::vector<size_t> buckets(bucketCount, 0);

  for (int64_t symbol : symbols) {
    ++buckets[static_cast<size_t>(symbol)];
  }

  return buckets;
}

// Fill `heads` with the start index of each bucket (for L-type induced sort).
void getBucketHeads(const std::vector<size_t> &buckets, std::vector<size_t> &heads) {
  size_t offset = 0;
  for (size_t i = 0; i < buckets.size(); ++i) {
    heads[i] = offset;
    offset += buckets[i];
  }
}

// Fill `tails` with the last index (inclusive) of each bucket (for S-type).
void getBucketTails(const std::vector<size_t> &buckets, std::vector<size_t> &tails) {
  size_t offset = 0;
  for (size_t i = 0; i < buckets.size(); ++i) {
    offset += buckets[i];
    tails[i] = offset - 1;
  }
}

// ── Induced sorting pass ───────────────────────────────────────────────────

// Forward pass: place all L-type suffixes using already-placed entries in sa.
void induceSortL(const std::vector<int64_t> &s, std::vector<size_t> &sa, const std::vector<bool> &isS,
                 const std::vector<size_t> &buckets, int64_t alphabetSize) {
  const size_t textSize = s.size();
  std::vector<size_t> heads(static_cast<size_t>(alphabetSize));
  getBucketHeads(buckets, heads);

  for (size_t i = 0; i < textSize; ++i) {
    const bool isEmptySlot = (sa[i] == SIZE_MAX);
    if (isEmptySlot) {
      continue;
    }

    const size_t suffixIndex = sa[i];
    if (suffixIndex == 0) {
      continue;
    }

    const size_t predecessor = suffixIndex - 1;
    if (isS[predecessor]) {
      continue;
    }

    const size_t bucket = static_cast<size_t>(s[predecessor]);
    sa[heads[bucket]++] = predecessor;
  }
}

// Backward pass: place all S-type suffixes.
void induceSortS(const std::vector<int64_t> &s, std::vector<size_t> &sa, const std::vector<bool> &isS,
                 const std::vector<size_t> &buckets, int64_t alphabetSize) {
  const size_t n = s.size();
  std::vector<size_t> tails(static_cast<size_t>(alphabetSize));
  getBucketTails(buckets, tails);

  for (size_t pos = n; pos-- > 0;) {
    const size_t suffix = sa[pos];
    if (suffix == SIZE_MAX || suffix == 0) {
      continue;
    }

    const size_t pred = suffix - 1;
    if (!isS[pred]) {
      continue;
    }

    const size_t bucket = static_cast<size_t>(s[pred]);
    sa[tails[bucket]--] = pred;
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// SuffixArray::sa_is  —  recursive core
// ─────────────────────────────────────────────────────────────────────────────

void SuffixArray::sa_is(const std::vector<int64_t> &s, std::vector<size_t> &sa,
                        int64_t alphabetSize) {
  constexpr size_t kEmptySlot = SIZE_MAX;
  const size_t n = s.size();

  std::vector<bool> isSType(n, false);
  isSType[n - 1] = true;
  for (size_t i = n - 1; i-- > 0;) {
    isSType[i] = (s[i] < s[i + 1]) || (s[i] == s[i + 1] && isSType[i + 1]);
  }

  auto isLMSPosition = [&](size_t i) -> bool {
    return i > 0 && i < n && isSType[i] && !isSType[i - 1];
  };

  const std::vector<size_t> buckets = getBuckets(s, alphabetSize);

  auto placeLMSIntoBuckets = [&](const std::vector<size_t> &orderedLMSPositions) {
        sa.assign(n, kEmptySlot);

        std::vector<size_t> tails(static_cast<size_t>(alphabetSize));
        getBucketTails(buckets, tails);

        for (size_t i = orderedLMSPositions.size() - 1; i > 0; --i) {
          const size_t pos = orderedLMSPositions[i];
          sa[tails[static_cast<size_t>(s[pos])]--] = pos;
        }

        sa[0] = n - 1;
      };

  std::vector<size_t> initialLMSPositions;
  initialLMSPositions.reserve(n / 2 + 1);
  for (size_t i = 0; i < n; ++i) {
    if (isLMSPosition(i)) {
      initialLMSPositions.push_back(i);
    }
  }

  placeLMSIntoBuckets(initialLMSPositions);
  induceSortL(s, sa, isSType, buckets, alphabetSize);
  induceSortS(s, sa, isSType, buckets, alphabetSize);

  std::vector<size_t> sortedLMSPositions;
  sortedLMSPositions.reserve(n / 2 + 1);
  for (size_t i = 0; i < n; ++i) {
    if (isLMSPosition(sa[i])) {
      sortedLMSPositions.push_back(sa[i]);
    }
  }

  if (sortedLMSPositions.empty()) {
    // No non-sentinel LMS positions; induced sorting already produced the SA.
    return;
  }

  std::vector<int64_t> rank(n, -1);
  int64_t currentRank = 0;
  rank[sortedLMSPositions[0]] = 0;

  for (size_t i = 1; i < sortedLMSPositions.size(); ++i) {
    const size_t previous = sortedLMSPositions[i - 1];
    const size_t current = sortedLMSPositions[i];

    bool differ = false;
    for (size_t d = 0;; ++d) { // loop length unknown
      const size_t prevIndex = previous + d;
      const size_t currIndex = current + d;

      if (prevIndex >= n || currIndex >= n) {
        differ = true;
        break;
      }

      const bool previousLMS = (d > 0) && isLMSPosition(prevIndex);
      const bool currentLMS = (d > 0) && isLMSPosition(currIndex);

      if (s[prevIndex] != s[currIndex] ||
          isSType[prevIndex] != isSType[currIndex] ||
          previousLMS != currentLMS) {
        differ = true;
        break;
      }

      if (d > 0 && (previousLMS || currentLMS)) {
        break;
      }
    }

    if (differ) {
      ++currentRank;
    }
    rank[current] = currentRank;
  }

  std::vector<size_t> lmsPositionsInTextOrder;
  lmsPositionsInTextOrder.reserve(sortedLMSPositions.size());
  for (size_t i = 0; i < n; ++i) {
    if (rank[i] >= 0) {
      lmsPositionsInTextOrder.push_back(i);
    }
  }

  std::vector<int64_t> reducedString;
  reducedString.reserve(lmsPositionsInTextOrder.size());
  for (size_t pos : lmsPositionsInTextOrder) {
    reducedString.push_back(rank[pos]);
  }

  const int64_t reducedAlphabet = currentRank + 1;
  std::vector<size_t> saReduced(reducedString.size());

  if (reducedAlphabet < static_cast<int64_t>(reducedString.size())) {
    sa_is(reducedString, saReduced, reducedAlphabet);
  } else {
    for (size_t i = 0; i < reducedString.size(); ++i) {
      saReduced[static_cast<size_t>(reducedString[i])] = i;
    }
  }

  std::vector<size_t> orderedLMSPositions(saReduced.size());
  for (size_t i = 0; i < saReduced.size(); ++i) {
    orderedLMSPositions[i] = lmsPositionsInTextOrder[saReduced[i]];
  }

  placeLMSIntoBuckets(orderedLMSPositions);
  induceSortL(s, sa, isSType, buckets, alphabetSize);
  induceSortS(s, sa, isSType, buckets, alphabetSize);
}

// ─────────────────────────────────────────────────────────────────────────────
// buildSuffixArray  —  driver
// ─────────────────────────────────────────────────────────────────────────────

void SuffixArray::buildSuffixArray() {
  if (_num == 0) {
    _sa.clear();
    return;
  }

  constexpr int64_t sentinel = 0;
  constexpr int64_t alphabetSize = 256;

  const size_t totalSize = _num + 1; // text + sentinel
  std::vector<int64_t> textWithSentinel(totalSize);

  for (size_t i = 0; i < _num; ++i) {
    textWithSentinel[i] = static_cast<unsigned char>(_data[i]);
  }
  textWithSentinel[_num] = sentinel;

  std::vector<size_t> suffixArray(totalSize);
  sa_is(textWithSentinel, suffixArray, alphabetSize);

  _sa.assign(suffixArray.begin() + 1, suffixArray.end());
}

// ─────────────────────────────────────────────────────────────────────────────
// Search  —  O(m log n) via two binary searches
// ─────────────────────────────────────────────────────────────────────────────

std::vector<SearchResult>
SuffixArray::search(const std::string &pattern) const {
  if (!_ready || pattern.empty()) {
    return {};
  }

  const size_t patternLength = pattern.size();
  const size_t firstMatch = lowerBound(pattern);
  const size_t pastLastMatch = upperBound(pattern);

  if (firstMatch >= pastLastMatch) {
    return {};
  }

  std::vector<SearchResult> results;
  results.reserve(pastLastMatch - firstMatch);

  for (size_t index = firstMatch; index < pastLastMatch; ++index) {
    results.emplace_back(SearchResult{_sa[index], patternLength});
  }

  std::sort(results.begin(), results.end(),
            [](const SearchResult &lhs, const SearchResult &rhs) {
              return lhs.offset < rhs.offset;
            });

  return results;
}

// ─────────────────────────────────────────────────────────────────────────────
// Binary-search helpers
// ─────────────────────────────────────────────────────────────────────────────

size_t SuffixArray::lowerBound(const std::string &pattern) const {
  const size_t patternLength = pattern.size();
  size_t low = 0;
  size_t high = _num;

  while (low < high) {
    const size_t mid = low + (high - low) / 2;
    const size_t suffixStart = _sa[mid];
    const size_t suffixLength = _num - suffixStart;
    const size_t compareLength = std::min(patternLength, suffixLength);

    const int cmp =
        pattern.compare(0, patternLength, _data + suffixStart, compareLength);

    if (cmp <= 0) {
      high = mid;
    } else {
      low = mid + 1;
    }
  }
  return low;
}

size_t SuffixArray::upperBound(const std::string &pattern) const {
  const size_t patternLength = pattern.size();
  size_t low = 0;
  size_t high = _num;

  while (low < high) {
    const size_t mid = low + (high - low) / 2;
    const size_t suffixPos = _sa[mid];
    const size_t suffixLength = _num - suffixPos;
    const size_t compareLength = std::min(patternLength, suffixLength);

    const int cmp =
        pattern.compare(0, patternLength, _data + suffixPos, compareLength);

    if (cmp < 0) {
      high = mid;
    } else {
      low = mid + 1;
    }
  }
  return low;
}
