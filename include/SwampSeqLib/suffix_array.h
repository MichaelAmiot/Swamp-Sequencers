#pragma once
#include "genome_mapper.h"
#include <cstddef>
#include <string>
#include <vector>

/*
  Result of a single pattern search.
 */
struct SearchResult {
  size_t offset; // Byte offset into the genome where the match starts.
  size_t length; // Length of the matched pattern (== query length).
};

/*
  Builds and queries a suffix array over a GenomeMapper view.

  Construction is O(n log n) time and O(n) extra space (the suffix array
  itself — the genome is never copied).  Queries are O(m log n).
 */
class SuffixArray {
public:
  /*
    Construct from an already-opened GenomeMapper.

    throws std::runtime_error if mapper.isValid() == false.
   */
  explicit SuffixArray(const GenomeMapper &mapper);
  SuffixArray(const std::string &text);

  /*
    Find all occurrences of @p pattern in the genome.

   */
  std::vector<SearchResult> search(const std::string &pattern) const;

  /*
   Return the number of suffixes (== genome size).
   */
  size_t size() const noexcept { return _n; }

  /*
   True when the suffix array was built successfully.
   */
  bool isReady() const noexcept { return _ready; }

private:
  // ──  helpers ────────────────────────────────────────────────────

  // Build the suffix array using prefix-doubling (Manber & Myers).
  void buildSuffixArray();

  /*
    Binary-search the SA for the leftmost suffix that begins with
    pattern.
   */
  size_t lowerBound(const std::string &pattern) const;

  /*
    Binary-search the SA for one-past the rightmost suffix that
     begins with pattern.
   */
  size_t upperBound(const std::string &pattern) const;

  // ── data members ────────────────────────────────────────────────────────

  const char *_data = nullptr; // Pointer into mapped genome memory.
  size_t _num = 0;             // Number of characters (genome length).
  std::vector<size_t> _sa;     // Suffix array: _sa[i] = start of i-th suffix.
  bool _ready = false;
  size_t _n = 0; // Count of found suffixes
};
