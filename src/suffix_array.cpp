#include "SwampSeqLib/suffix_array.h"
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <limits>

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
    _owned = text;
    _data = _owned.c_str();
    _num = _owned.size();

    buildSuffixArray();
    _ready = true;
}

// ─────────────────────────────────────────────────────────────────────────────
// SA-IS Helpers
// ─────────────────────────────────────────────────────────────────────────────

// ── Helpers (file-scope) ─────────────────────────────────

// Return bucket sizes: buckets[c] = number of occurrences of symbol c.
std::vector<size_t> getBuckets(const std::vector<uint32_t> &symbols, uint32_t alphabetSize) {
    // Safety check: Ensure bucket size covers the requested alphabet size
    std::vector<size_t> buckets(alphabetSize, 0);

    for (uint32_t symbol : symbols) {
        buckets[static_cast<size_t>(symbol)]++;
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

void induceSortL(const std::vector<uint32_t> &s, std::vector<size_t> &sa, const std::vector<bool> &isSType,
                 const std::vector<size_t> &buckets, uint32_t alphabetSize) {
    const size_t textSize = s.size();
    std::vector<size_t> heads(alphabetSize);
    getBucketHeads(buckets, heads);

    for (size_t i = 0; i < textSize; ++i) {
        if (sa[i] == SIZE_MAX || sa[i] == 0) continue;

        // Place the predecessor at the current head of its bucket.
        const size_t predecessor = sa[i] - 1;
        if (!isSType[predecessor]) {
            sa[heads[s[predecessor]]++] = predecessor;
        }
    }
}

void induceSortS(const std::vector<uint32_t> &s, std::vector<size_t> &sa, const std::vector<bool> &isSType,
                 const std::vector<size_t> &buckets, uint32_t alphabetSize) {
    const size_t n = s.size();
    std::vector<size_t> tails(alphabetSize);
    getBucketTails(buckets, tails);

    for (size_t i = n; i-- > 0;) {
        if (sa[i] == SIZE_MAX || sa[i] == 0) continue;

        const size_t predecessor = sa[i] - 1;
        if (isSType[predecessor]) {
            sa[tails[s[predecessor]]--] = predecessor;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SuffixArray::sa_is
// ─────────────────────────────────────────────────────────────────────────────

void SuffixArray::sa_is(const std::vector<uint32_t> &s, std::vector<size_t> &sa, uint32_t alphabetSize) {
    const size_t kUnusedSlot = SIZE_MAX;
    const uint32_t kEmptyRank = std::numeric_limits<uint32_t>::max();
    const size_t n = s.size();

    std::vector<bool> isSType(n, false);
    isSType[n - 1] = true;
    for (size_t i = n - 1; i-- > 0;) {
        isSType[i] = (s[i] < s[i + 1]) || (s[i] == s[i + 1] && isSType[i + 1]);
    }

    auto isLMSPosition = [&](size_t i) -> bool {
        return i > 0 && isSType[i] && !isSType[i - 1];
    };

    const std::vector<size_t> buckets = getBuckets(s, alphabetSize);

    auto placeLMSIntoBuckets = [&](const std::vector<size_t> &orderedLMSPositions) {
        sa.assign(n, kUnusedSlot);
        std::vector<size_t> tails(alphabetSize);
        getBucketTails(buckets, tails);

        for (size_t i = orderedLMSPositions.size(); i-- > 0;) {
            const size_t pos = orderedLMSPositions[i];
            sa[tails[s[pos]]--] = pos;
        }
        sa[0] = n - 1; // Sentinel
    };

    std::vector<size_t> initialLMSPositions;
    for (size_t i = 0; i < n; ++i) {
        if (isLMSPosition(i)) initialLMSPositions.push_back(i);
    }

    placeLMSIntoBuckets(initialLMSPositions);
    induceSortL(s, sa, isSType, buckets, alphabetSize);
    induceSortS(s, sa, isSType, buckets, alphabetSize);

    std::vector<size_t> sortedLMSPositions;
    for (size_t i = 0; i < n; ++i) {
        if (isLMSPosition(sa[i])) sortedLMSPositions.push_back(sa[i]);
    }

    if (sortedLMSPositions.empty()) return;

    // FIX: Using a sentinel that doesn't wrap to 0 or overflow logic
    std::vector<uint32_t> rank(n, kEmptyRank);
    uint32_t currentRank = 0;
    rank[sortedLMSPositions[0]] = 0;

    for (size_t i = 1; i < sortedLMSPositions.size(); ++i) {
        const size_t previous = sortedLMSPositions[i - 1];
        const size_t current = sortedLMSPositions[i];
        bool differ = false;

        for (size_t d = 0;; ++d) {
            const size_t pIdx = previous + d;
            const size_t cIdx = current + d;

            const bool pLMS = (d > 0) && isLMSPosition(pIdx);
            const bool cLMS = (d > 0) && isLMSPosition(cIdx);

            if (s[pIdx] != s[cIdx] || isSType[pIdx] != isSType[cIdx] || pLMS != cLMS) {
                differ = true;
                break;
            }
            if (d > 0 && (pLMS || cLMS)) break;
        }

        if (differ) ++currentRank;
        rank[current] = currentRank;
    }

    std::vector<size_t> lmsPositionsInTextOrder;
    std::vector<uint32_t> reducedString;
    for (size_t i = 0; i < n; ++i) {
        // FIX: Explicitly check against the kEmptyRank sentinel
        if (rank[i] != kEmptyRank) {
            lmsPositionsInTextOrder.push_back(i);
            reducedString.push_back(rank[i]);
        }
    }

    const uint32_t reducedAlphabetSize = currentRank + 1;
    std::vector<size_t> saReduced(reducedString.size());

    if (reducedAlphabetSize < static_cast<uint32_t>(reducedString.size())) {
        sa_is(reducedString, saReduced, reducedAlphabetSize);
    } else {
        for (size_t i = 0; i < reducedString.size(); ++i) {
            saReduced[reducedString[i]] = i;
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
// SuffixArray Core Logic
// ─────────────────────────────────────────────────────────────────────────────

void SuffixArray::buildSuffixArray() {
    if (_num == 0) {
        _sa.clear();
        return;
    }

    constexpr uint32_t alphabetSize = 256;
    const size_t totalSize = _num + 1;
    std::vector<uint32_t> textWithSentinel(totalSize);

    for (size_t i = 0; i < _num; ++i) {
        textWithSentinel[i] = static_cast<unsigned char>(_data[i]);
    }
    textWithSentinel[_num] = 0; // Sentinel

    std::vector<size_t> suffixArray(totalSize);
    sa_is(textWithSentinel, suffixArray, alphabetSize);

    _sa.assign(suffixArray.begin() + 1, suffixArray.end());
}

std::vector<SearchResult> SuffixArray::search(const std::string &pattern) const {
    if (!_ready || pattern.empty()) return {};

    const size_t firstMatch = lowerBound(pattern);
    const size_t pastLastMatch = upperBound(pattern);

    if (firstMatch >= pastLastMatch) return {};

    std::vector<SearchResult> results;
    results.reserve(pastLastMatch - firstMatch);

    for (size_t index = firstMatch; index < pastLastMatch; ++index) {
        results.emplace_back(SearchResult{_sa[index], (uint32_t)pattern.size()});
    }

    std::sort(results.begin(), results.end(),
              [](const SearchResult &lhs, const SearchResult &rhs) {
                  return lhs.offset < rhs.offset;
              });

    return results;
}

size_t SuffixArray::size() const noexcept {
    return _num;
}

bool SuffixArray::ready() const noexcept {
    return _ready;
}

const std::vector<size_t> & SuffixArray::sa() const noexcept {
    return _sa;
}

size_t SuffixArray::lowerBound(const std::string &pattern) const {
    const size_t patternLength = pattern.size();
    size_t low = 0, high = _num;

    while (low < high) {
        const size_t mid = low + (high - low) / 2;
        const size_t suffixStart = _sa[mid];
        const size_t suffixLength = _num - suffixStart;
        const size_t compareLength = std::min(patternLength, suffixLength);

        const int cmp = pattern.compare(0, patternLength, _data + suffixStart, compareLength);

        if (cmp <= 0) high = mid;
        else low = mid + 1;
    }
    return low;
}

size_t SuffixArray::upperBound(const std::string &pattern) const {
    const size_t patternLength = pattern.size();
    size_t low = 0, high = _num;

    while (low < high) {
        const size_t mid = low + (high - low) / 2;
        const size_t suffixPos = _sa[mid];
        const size_t suffixLength = _num - suffixPos;
        const size_t compareLength = std::min(patternLength, suffixLength);

        const int cmp = pattern.compare(0, patternLength, _data + suffixPos, compareLength);

        if (cmp < 0) high = mid;
        else low = mid + 1;
    }
    return low;
}