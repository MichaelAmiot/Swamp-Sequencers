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
//     the original string, yielding the final Suffix Array.
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
std::vector<size_t> getBuckets(const std::vector<uint32_t> &symbols, uint32_t alphabetSize) {
    const size_t bucketCount = alphabetSize;
    std::vector<size_t> buckets(bucketCount, 0);

    for (uint32_t symbol: symbols) {
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
void induceSortL(const std::vector<uint32_t> &s, std::vector<size_t> &sa, const std::vector<bool> &isSType,
                 const std::vector<size_t> &buckets, uint32_t alphabetSize) {
    const size_t textSize = s.size();
    std::vector<size_t> heads(alphabetSize);
    getBucketHeads(buckets, heads);

    for (size_t i = 0; i < textSize; ++i) {
        const size_t suffixIndex = sa[i];

        // Skip empty slots in the suffix array.
        if (suffixIndex == SIZE_MAX) {
            continue;
        }

        // The sentinel suffix has no predecessor.
        if (suffixIndex == 0) {
            continue;
        }

        const size_t predecessor = suffixIndex - 1;

        // Only L-type predecessors are induced in this pass.
        if (isSType[predecessor]) {
            continue;
        }

        // Place the predecessor at the current head of its bucket.
        const size_t bucket = s[predecessor];
        sa[heads[bucket]++] = predecessor;
    }
}

// Backward pass: place all S-type suffixes.
void induceSortS(const std::vector<uint32_t> &s, std::vector<size_t> &sa, const std::vector<bool> &is_S,
                 const std::vector<size_t> &buckets, uint32_t alphabetSize) {
    const size_t n = s.size();
    std::vector<size_t> tails(alphabetSize);
    getBucketTails(buckets, tails);

    for (size_t pos = n; pos-- > 0;) {
        const size_t suffix = sa[pos];
        if (suffix == SIZE_MAX || suffix == 0) {
            continue;
        }

        const size_t pred = suffix - 1;
        if (!is_S[pred]) {
            continue;
        }

        const size_t bucket = s[pred];
        sa[tails[bucket]--] = pred;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SuffixArray::sa_is  —  recursive 
// ─────────────────────────────────────────────────────────────────────────────

void SuffixArray::sa_is(const std::vector<uint32_t> &s, std::vector<size_t> &sa, uint32_t alphabetSize) {
    const size_t kUnusedSlot = SIZE_MAX;
    const size_t n = s.size();

    // Classify every position as S-type or L-type.
    // The last position is always S-type because it is the sentinel.
    std::vector isSType(n, false);
    isSType[n - 1] = true;
    for (size_t i = n - 1; i-- > 0;) {
        isSType[i] = (s[i] < s[i + 1]) || (s[i] == s[i + 1] && isSType[i + 1]);
    }

    // LMS positions are S-type positions whose predecessor is L-type.
    auto isLMSPosition = [&](size_t i) -> bool {
        return i > 0 && i < n && isSType[i] && !isSType[i - 1];
    };

    const std::vector<size_t> buckets = getBuckets(s, alphabetSize);

    // Place LMS positions into the tail of each bucket.
    auto placeLMSIntoBuckets = [&](const std::vector<size_t> &orderedLMSPositions) {
        sa.assign(n, kUnusedSlot);

        std::vector<size_t> tails(alphabetSize);
        getBucketTails(buckets, tails);

        // Insert from right to left so relative bucket order is preserved.
        for (size_t i = orderedLMSPositions.size(); i-- > 0;) {
            const size_t pos = orderedLMSPositions[i];
            const auto bucketIndex = static_cast<size_t>(s[pos]);
            sa[tails[bucketIndex]--] = pos;
        }

        // Sentinel suffix is always the smallest suffix.
        sa[0] = n - 1;
    };

    // Collect LMS positions in text order for the initial induced sort.
    std::vector<size_t> initialLMSPositions;
    initialLMSPositions.reserve(n / 2 + 1);
    for (size_t i = 0; i < n; ++i) {
        if (isLMSPosition(i)) {
            initialLMSPositions.push_back(i);
        }
    }

    // First induced sort pass using the LMS suffixes in text order.
    placeLMSIntoBuckets(initialLMSPositions);
    induceSortL(s, sa, isSType, buckets, alphabetSize);
    induceSortS(s, sa, isSType, buckets, alphabetSize);

    // Extract LMS suffixes in their induced-sorted order.
    std::vector<size_t> sortedLMSPositions;
    sortedLMSPositions.reserve(n / 2 + 1);
    for (size_t i = 0; i < n; ++i) {
        if (isLMSPosition(sa[i])) {
            sortedLMSPositions.push_back(sa[i]);
        }
    }

    // If there are no non-sentinel LMS positions, the induced sort is complete.
    if (sortedLMSPositions.empty()) {
        return;
    }

    // Assign ranks to LMS substrings.
    // Equal substrings receive the same rank, different ones increment the rank.
    std::vector<uint32_t> rank(n, -1);
    uint32_t currentRank = 0;
    rank[sortedLMSPositions[0]] = 0;

    // Compare adjacent LMS substrings to decide whether their ranks differ.
    for (size_t i = 1; i < sortedLMSPositions.size(); ++i) {
        const size_t previous = sortedLMSPositions[i - 1];
        const size_t current = sortedLMSPositions[i];

        bool differ = false;

        // Compare the two LMS substrings character by character.
        // The loop ends when either substring reaches the next LMS boundary
        // or when a difference is found.
        for (size_t d = 0;; ++d) {
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

            // Stop once both substrings reach their LMS boundary.
            if (d > 0 && (previousLMS || currentLMS)) {
                break;
            }
        }

        if (differ) {
            ++currentRank;
        }
        rank[current] = currentRank;
    }

    // Build the reduced string from LMS ranks in text order.
    std::vector<size_t> lmsPositionsInTextOrder;
    lmsPositionsInTextOrder.reserve(sortedLMSPositions.size());
    for (size_t i = 0; i < n; ++i) {
        if (rank[i] >= 0) {
            lmsPositionsInTextOrder.push_back(i);
        }
    }

    std::vector<uint32_t> reducedString;
    reducedString.reserve(lmsPositionsInTextOrder.size());
    for (size_t pos: lmsPositionsInTextOrder) {
        reducedString.push_back(rank[pos]);
    }

    // Recurse only if the reduced alphabet is still smaller than the string.
    const uint32_t reducedAlphabet = currentRank + 1;
    std::vector<size_t> saReduced(reducedString.size());

    if (reducedAlphabet < static_cast<uint32_t>(reducedString.size())) {
        sa_is(reducedString, saReduced, reducedAlphabet);
    } else {
        // Directly map each rank to its position when the reduced problem is trivial.
        for (size_t i = 0; i < reducedString.size(); ++i) {
            saReduced[static_cast<size_t>(reducedString[i])] = i;
        }
    }

    // Translate reduced-string suffix order back to original LMS positions.
    std::vector<size_t> orderedLMSPositions(saReduced.size());
    for (size_t i = 0; i < saReduced.size(); ++i) {
        orderedLMSPositions[i] = lmsPositionsInTextOrder[saReduced[i]];
    }

    // Final induced sort using the exact LMS order recovered from the reduced problem.
    placeLMSIntoBuckets(orderedLMSPositions);
    induceSortL(s, sa, isSType, buckets, alphabetSize);
    induceSortS(s, sa, isSType, buckets, alphabetSize);
}

// ─────────────────────────────────────────────────────────────────────────────
// buildSuffixArray
// ─────────────────────────────────────────────────────────────────────────────

void SuffixArray::buildSuffixArray() {
    if (_num == 0) {
        _sa.clear();
        return;
    }

    constexpr uint32_t sentinel = 0;
    constexpr uint32_t alphabetSize = 256;

    // Encode the input text as integer symbols and append a unique sentinel.
    const size_t totalSize = _num + 1; // text + sentinel
    std::vector<uint32_t> textWithSentinel(totalSize);

    for (size_t i = 0; i < _num; ++i) {
        textWithSentinel[i] = static_cast<unsigned char>(_data[i]);
    }
    textWithSentinel[_num] = sentinel;

    // Build the suffix array over the encoded text.
    std::vector<size_t> suffixArray(totalSize);
    sa_is(textWithSentinel, suffixArray, alphabetSize);

    // Remove the sentinel suffix entry from the final result.
    _sa.assign(suffixArray.begin() + 1, suffixArray.end());
}

// ─────────────────────────────────────────────────────────────────────────────
// Search  —  O(m log n) via two binary searches
// ─────────────────────────────────────────────────────────────────────────────

std::vector<SearchResult>
SuffixArray::search(const std::string &pattern) const {
    // Return early when the suffix array is not ready or the pattern is empty.
    if (!_ready || pattern.empty()) {
        return {};
    }

    // Find the inclusive/exclusive range of matching suffixes in the array.
    const size_t patternLength = pattern.size();
    const size_t firstMatch = lowerBound(pattern);
    const size_t pastLastMatch = upperBound(pattern);

    if (firstMatch >= pastLastMatch) {
        return {};
    }

    // Convert the matching suffix-array slice into SearchResult entries.
    std::vector<SearchResult> results;
    results.reserve(pastLastMatch - firstMatch);

    for (size_t index = firstMatch; index < pastLastMatch; ++index) {
        results.emplace_back(SearchResult{_sa[index], patternLength});
    }

    // Keep results sorted by original text offset for predictable output.
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
    // Search for the first suffix that is not lexicographically smaller than pattern.
    const size_t patternLength = pattern.size();
    size_t low = 0;
    size_t high = _num;

    while (low < high) {
        // Standard binary-search midpoint.
        const size_t mid = low + (high - low) / 2;

        // Compare the pattern with the suffix starting at _sa[mid].
        const size_t suffixStart = _sa[mid];
        const size_t suffixLength = _num - suffixStart;
        const size_t compareLength = std::min(patternLength, suffixLength);
        const int cmp = pattern.compare(0, patternLength, _data + suffixStart, compareLength);

        // If suffix >= pattern, keep the left half; otherwise move right.
        if (cmp <= 0) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }

    return low;
}

size_t SuffixArray::upperBound(const std::string &pattern) const {
    // Search for the first suffix that is lexicographically greater than the pattern.
    const size_t patternLength = pattern.size();
    size_t low = 0;
    size_t high = _num;

    while (low < high) {
        // Standard binary-search midpoint calculation.
        const size_t mid = low + (high - low) / 2;

        // Compare the pattern against the suffix starting at _sa[mid].
        const size_t suffixPos = _sa[mid];
        const size_t suffixLength = _num - suffixPos;
        const size_t compareLength = std::min(patternLength, suffixLength);
        const int cmp = pattern.compare(0, patternLength, _data + suffixPos, compareLength);

        // Move left when the suffix is already greater than the pattern;
        // otherwise, continue searching in the right half.
        if (cmp < 0) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }

    return low;
}