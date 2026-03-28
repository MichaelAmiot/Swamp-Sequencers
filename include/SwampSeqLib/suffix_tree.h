#pragma once
#include "SwampSeqLib/genome_mapper.h"
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// STSearchResult
//   Returned by SuffixTree::search().
//   offset – 0-based position of the match in the original text.
//   length – length of the matched pattern (equals pattern.size()).
// ─────────────────────────────────────────────────────────────────────────────
struct STSearchResult {
  int64_t offset;
  int64_t length;
  bool operator==(const STSearchResult &other) const;
};

// ─────────────────────────────────────────────────────────────────────────────
// SuffixTree
//
// Builds a compressed trie (Patricia trie) of all suffixes of a text using
// Ukkonen's online algorithm in O(n) time and O(n) space.
//
// Every internal node and leaf is represented by a Node struct.  Edges are
// labelled by half-open intervals [start, end) into the original text, so no
// substring copies are made during construction.
//
// Leaf end pointers all point to a single shared "global end" int64_t that is
// incremented with each character — this is the core of Ukkonen's O(n)
// guarantee (open-ended leaves never need updating individually).
//
// Lifetime rules mirror SuffixArray:
//   • GenomeMapper constructor — raw pointer stored; mapper must outlive this.
//   • std::string constructor  — text copied into _owned; self-managed.
//
// Search runs in O(m) where m = |pattern|.
// ─────────────────────────────────────────────────────────────────────────────
class SuffixTree {
public:
  // ── Construction ──────────────────────────────────────────────────────────

  // Build from a GenomeMapper. Mapper must outlive this object.
  explicit SuffixTree(const GenomeMapper &mapper);

  // Build from an arbitrary string (text copied internally).
  explicit SuffixTree(const std::string &text);

  // Non-copyable; movable.
  SuffixTree(const SuffixTree &)            = delete;
  SuffixTree &operator=(const SuffixTree &) = delete;
  SuffixTree(SuffixTree &&)                 = default;
  SuffixTree &operator=(SuffixTree &&)      = default;

  ~SuffixTree() = default;

  // ── Query ─────────────────────────────────────────────────────────────────

  // Return all positions at which pattern occurs in the text, sorted by
  // ascending offset. Returns empty vector if pattern is absent or the tree
  // has not been built successfully.
  std::vector<STSearchResult> search(const std::string &pattern) const;

  // Number of characters in the indexed text.
  int64_t size() const noexcept { return _num; }

  // True once the suffix tree has been successfully built.
  bool ready() const noexcept { return _ready; }

private:
  // ── Node ──────────────────────────────────────────────────────────────────
  //
  // Each node represents a state in the implicit suffix tree.
  //
  //  start      – index into _text where this edge label begins.
  //  end        – pointer to the (exclusive) end of the edge label.
  //               For leaves this points to _globalEnd (shared, open).
  //               For internal nodes this points to _nodeEnds[nodeIndex].
  //  suffixLink – Ukkonen suffix link; points to the longest proper suffix
  //               state reachable from this node (-1 = unset).
  //  suffixIndex – for leaves: which suffix starts here (-1 for internal).
  //  children   – map from first character of child edge → child node index.
  //               Using unordered_map keeps child lookup O(1) average.

  struct Node {
    int64_t  start      = -1;
    int64_t *end        = nullptr;   // points into _nodeEnds or _globalEnd
    int64_t  suffixLink = -1;
    int64_t  suffixIndex= -1;        // -1 for internal nodes
    std::unordered_map<int64_t, int64_t> children;

    // Edge length. Leaves use *end which equals _globalEnd at query time.
    int64_t edgeLength() const { return *end - start; }
  };

  // ── Internal data ─────────────────────────────────────────────────────────

  const char  *_data  = nullptr;   // raw text pointer (not owned for mapper)
  std::string  _owned;             // owned copy when built from std::string
  int64_t      _num   = 0;         // text length (excludes sentinel)
  bool         _ready = false;

  // Node pool — avoids pointer invalidation on reallocation.
  std::vector<Node>    _nodes;
  // Per-internal-node end values (leaves share _globalEnd).
  std::vector<int64_t> _nodeEnds;
  // The shared open end for all leaves; incremented each phase.
  int64_t              _globalEnd = 0;

  // Ukkonen active point
  int64_t _activeNode   = 0;   // index into _nodes
  int64_t _activeEdge   = -1;  // first character of the active edge (as index)
  int64_t _activeLength = 0;
  int64_t _remaining    = 0;   // suffixes yet to be inserted

  // Root is always node 0; -1 sentinel node for suffix links from root.
  static constexpr int64_t ROOT     =  0;
  static constexpr int64_t NO_NODE  = -1;

  // ── Ukkonen construction ──────────────────────────────────────────────────

  void buildSuffixTree();

  // Allocate a new leaf node for the edge starting at `start`.
  int64_t newLeaf(int64_t start);

  // Allocate a new internal node for the edge interval [start, end).
  int64_t newInternal(int64_t start, int64_t end);

  // Extend the tree by one character at position `pos`.
  void extendTree(int64_t pos);

  // ── Post-construction annotation ─────────────────────────────────────────

  // DFS to stamp suffix indices onto every leaf.
  void annotateSuffixIndices(int64_t nodeIdx, int64_t labelHeight);

  // ── Search helpers ────────────────────────────────────────────────────────

  // Collect all leaf suffixIndex values in the subtree rooted at nodeIdx.
  void collectLeaves(int64_t nodeIdx, std::vector<STSearchResult> &out,
                     int64_t patternLength) const;
};
