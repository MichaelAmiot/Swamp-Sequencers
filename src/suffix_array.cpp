#include "SwampSeqLib/suffix_array.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string_view>

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

// ─────────────────────────────────────────────────────────────────────────────
// Suffix Array Construction  —  O(n log² n) prefix-doubling (Manber & Myers)
// Simple, cache-friendly, no external dependencies.
// For genomes in the hundreds-of-MB range the O(n log²n) constant is
//  acceptable.  SA-IS (O(n)) can be substituted later if needed.

// Memory layout:
//   _sa[i]  = starting index of the i-th lexicographically smallest suffix.
// ─────────────────────────────────────────────────────────────────────────────

void SuffixArray::buildSuffixArray() {
  if (_num == 0) {
    _sa.clear();
    return;
  }

  // Initialise SA with identity permutation [0, 1, 2, …, n-1].
  _sa.resize(_num);
  std::iota(_sa.begin(), _sa.end(), size_t{0});

  // ranks[i] = current sort rank of suffix starting at i.
  std::vector<int> ranks(_num), nextRanks(_num);

  // Seed ranks from the first character.
  for (size_t i = 0; i < _num; ++i) {
    ranks[i] = static_cast<unsigned char>(_data[i]);
  }

  // Double the comparison window until the whole SA is sorted.
  for (size_t window = 1; window < _num; window <<= 1) {
    auto compareByIndex = [&](size_t leftIndex, size_t rightIndex) {
      if (ranks[leftIndex] != ranks[rightIndex]) {
        return ranks[leftIndex] < ranks[rightIndex];
      }

      const int rankA =
          (leftIndex + window < _num) ? ranks[leftIndex + window] : -1;
      const int rankB =
          (rightIndex + window < _num) ? ranks[rightIndex + window] : -1;
      return rankA < rankB;
    };
    std::sort(_sa.begin(), _sa.end(), compareByIndex);

    // Re-assign ranks based on the freshly sorted order.
    nextRanks[_sa[0]] = 0;
    for (size_t i = 1; i < _num; ++i) {
      nextRanks[_sa[i]] =
          nextRanks[_sa[i - 1]] + (compareByIndex(_sa[i - 1], _sa[i]) ? 1 : 0);
    }

    ranks = nextRanks;

    // Early exit: all ranks are unique → array is fully sorted.
    if (ranks[_sa[_n - 1]] == static_cast<int>(_num) - 1) {
      break;
    }
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Search  —  O(m log n) via two binary searches
// ─────────────────────────────────────────────────────────────────────────────

std::vector<SearchResult>
SuffixArray::search(const std::string &pattern) const {
  if (!_ready || pattern.empty())
    return {};

  const size_t patternLength = pattern.size();
  const size_t low = lowerBound(pattern);
  const size_t high = upperBound(pattern);

  if (low >= high)
    return {};

  std::vector<SearchResult> results;
  results.reserve(high - low);

  for (size_t i = low; i < high; ++i) {
    results.push_back({_sa[i], patternLength});
  }

  std::sort(results.begin(), results.end(),
            [](const SearchResult &left, const SearchResult &right) {
              return left.offset < right.offset;
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
