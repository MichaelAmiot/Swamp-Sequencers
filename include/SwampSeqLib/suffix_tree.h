#pragma once
#include "SwampSeqLib/genome_mapper.h"
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

struct STSearchResult {
    int64_t offset;
    int64_t length;
    bool operator==(const STSearchResult &other) const;
};

class SuffixTree {
public:
    explicit SuffixTree(const GenomeMapper &mapper);
    explicit SuffixTree(const std::string &text);

    SuffixTree(const SuffixTree &) = delete;
    SuffixTree &operator=(const SuffixTree &) = delete;
    SuffixTree(SuffixTree &&) = default;
    SuffixTree &operator=(SuffixTree &&) = default;
    ~SuffixTree() = default;

    std::vector<STSearchResult> search(const std::string &pattern) const;

    int64_t size()  const noexcept { return _textLength; }
    bool ready() const noexcept { return _isReady; }

private:

    static constexpr int kAlphabetSize = 6;
    static constexpr int64_t kInvalidIndex = std::numeric_limits<int64_t>::max();
    static constexpr int64_t kLeafEnd = kInvalidIndex;

    static int charToSlot(unsigned char c) noexcept;

    struct Node {
        int64_t start = 0;
        int64_t suffixIndex = kInvalidIndex;
        int64_t suffixLink  = 0;
        int64_t endIdx = kInvalidIndex;
        std::array<int64_t, kAlphabetSize> children{};

        int64_t edgeLength(int64_t currentGlobalEnd, const std::vector<int64_t> &nodeEnds) const noexcept;
        int64_t edgeEnd(int64_t currentGlobalEnd, const std::vector<int64_t> &nodeEnds) const noexcept;
    };

    const char  *_textData = nullptr;
    std::string _owned;
    int64_t _textLength = 0;
    bool _isReady = false;
    std::vector<Node>_nodes;
    std::vector<int64_t>_nodeEnds;
    int64_t _currentGlobalEnd = 0;

    int64_t _activeNode   = 0;
    int _activeSlot   = -1;
    int64_t _activeLength = 0;
    int64_t _remaining    = 0;

    static constexpr int64_t kRoot = 0;

    void buildSuffixTree();
    void initializeBuildState();
    int64_t newLeaf (int64_t start, int64_t suffixIndex);
    int64_t newInternal(int64_t start, int64_t end);
    void extendTree (int64_t pos);
   // void annotateSuffixIndices();

    void collectLeaves(int64_t subtreeRoot, std::vector<STSearchResult> &out, int64_t patternLength) const;
};