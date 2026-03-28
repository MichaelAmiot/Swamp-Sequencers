#pragma once
#include "SwampSeqLib/genome_mapper.h"
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// STSearchResult
//   offset – 0-based position of the match in the original text.
//   length – length of the matched pattern.
// ─────────────────────────────────────────────────────────────────────────────
struct STSearchResult {
    uint32_t offset;
    uint32_t length;

    bool operator==(const STSearchResult &other) const;
};

// ─────────────────────────────────────────────────────────────────────────────
// SuffixTree

// Ukkonen's online O(n) suffix tree, optimised for genomic alphabets.

// Hard cap: 500 Mbp. Above that, peak RAM exceeds ~20 GB which is
// impractical for a suffix tree; use SuffixArray instead.

// Alphabet mapping (ALPHA = 6):
//   A→0  C→1  G→2  T→3  N→4  sentinel→5
//   Any other byte maps to slot 4 (N) as a safe fallback.

// Lifetime:
//   GenomeMapper constructor — raw pointer stored; mapper must outlive this.
//   std::string constructor  — text copied into _owned; self-managed.

// Search: O(m), m = pattern length.
// ─────────────────────────────────────────────────────────────────────────────
class SuffixTree {
public:
    // ── Construction ──────────────────────────────────────────────────────────
    explicit SuffixTree(const GenomeMapper &mapper);

    explicit SuffixTree(const std::string &text);

    SuffixTree(const SuffixTree &) = delete;

    SuffixTree &operator=(const SuffixTree &) = delete;

    SuffixTree(SuffixTree &&) = default;

    SuffixTree &operator=(SuffixTree &&) = default;

    ~SuffixTree() = default;

    // ── Query ─────────────────────────────────────────────────────────────────
    // All positions where pattern occurs, sorted ascending. Empty if absent.
    std::vector<STSearchResult> search(const std::string &pattern) const;

    uint32_t size() const noexcept { return _textLength; }
    bool ready() const noexcept { return _isReady; }

private:
    // ── Alphabet ───────────────────────────────────────────────────────────────
    static constexpr int kAlphabetSize = 6; // A C G T N sentinel
    static constexpr uint32_t kInvalidIndex = std::numeric_limits<uint32_t>::max();
    static constexpr uint32_t kLeafEnd = kInvalidIndex; // endIdx value for leaves

    // Map a text byte to a child-array slot [0, kAlphabetSize).
    static int charToSlot(unsigned char c) noexcept;

    // ── Node (40 bytes) ───────────────────────────────────────────────────────
    struct Node {
        uint32_t start = 0;
        uint32_t suffixIndex = kInvalidIndex; // kInvalidIndex for internal nodes
        uint32_t suffixLink = 0; // default to root
        uint32_t endIdx = kInvalidIndex; // kLeafEnd for leaves; else index into _nodeEnds
        std::array<uint32_t, kAlphabetSize> children{};

        uint32_t edgeLength(uint32_t currentGlobalEnd, const std::vector<uint32_t> &nodeEnds) const noexcept;
        uint32_t edgeEnd(uint32_t currentGlobalEnd, const std::vector<uint32_t> &nodeEnds) const noexcept;

    };

    // ── Internal data ─────────────────────────────────────────────────────────
    const char *_textData = nullptr;
    std::string _owned;
    uint32_t _textLength = 0;
    bool _isReady = false;
    std::vector<Node> _nodes;
    std::vector<uint32_t> _nodeEnds; // end values for internal nodes
    uint32_t _currentGlobalEnd = 0;

    // Ukkonen active point
    uint32_t _activeNode = 0;
    int _activeSlot = -1; // child slot of the active edge
    uint32_t _activeLength = 0;
    uint32_t _remaining = 0;

    static constexpr uint32_t kRoot = 0;

    // ── Build ─────────────────────────────────────────────────────────────────
    void buildSuffixTree();

    void initializeBuildState();

    uint32_t newLeaf(uint32_t start);

    uint32_t newInternal(uint32_t start, uint32_t end);

    void extendTree(uint32_t pos);

    void annotateSuffixIndices(); // iterative DFS

    // ── Search ────────────────────────────────────────────────────────────────
    void collectLeaves(uint32_t subtreeRoot, std::vector<STSearchResult> &out, uint32_t patternLength) const; // iterative DFS
};