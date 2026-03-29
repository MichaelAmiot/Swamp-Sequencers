#include "SwampSeqLib/suffix_tree.h"
#include <algorithm>
#include <cstdint>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────
// STSearchResult
// ─────────────────────────────────────────────────────────────────────────────

bool STSearchResult::operator==(const STSearchResult &other) const {
    return (this->offset == other.offset) && (this->length == other.length);
}

// ─────────────────────────────────────────────────────────────────────────────
// charToSlot
// ─────────────────────────────────────────────────────────────────────────────

int SuffixTree::charToSlot(unsigned char c) noexcept {
    switch (c) {
        case 'A':
        case 'a': return 0;
        case 'C':
        case 'c': return 1;
        case 'G':
        case 'g': return 2;
        case 'T':
        case 't': return 3;
        case 255: return 5; // sentinel
        default: return 4; // N and anything else
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Construction
// ─────────────────────────────────────────────────────────────────────────────

SuffixTree::SuffixTree(const GenomeMapper &mapper) {
    if (!mapper.isValid())
        throw std::runtime_error("SuffixTree: GenomeMapper is not valid.");

    constexpr int64_t MAX_BP = 2'000'000'000LL;
    if (static_cast<int64_t>(mapper.size()) > MAX_BP)
        throw std::runtime_error(
            "SuffixTree: genome exceeds 2 Gbp hard cap. "
            "Use SuffixArray for large genomes.");

    _textData = mapper.data();
    _textLength = static_cast<int64_t>(mapper.size());

    buildSuffixTree();
    _isReady = true;
}

SuffixTree::SuffixTree(const std::string &text) {
    constexpr int64_t MAX_BP = 2'000'000'000LL;
    if (static_cast<int64_t>(text.size()) > MAX_BP)
        throw std::runtime_error(
            "SuffixTree: text exceeds 2 Gbp hard cap. "
            "Use SuffixArray for large genomes.");

    _owned = text;
    _textData = _owned.c_str();
    _textLength = static_cast<int64_t>(_owned.size());

    buildSuffixTree();
    _isReady = true;
}

// ─────────────────────────────────────────────────────────────────────────────
// Node helpers
// ─────────────────────────────────────────────────────────────────────────────

int64_t SuffixTree::Node::edgeLength(int64_t currentGlobalEnd, const std::vector<int64_t> &nodeEnds) const noexcept {
    const int64_t end = (endIdx == kLeafEnd) ? currentGlobalEnd : nodeEnds[endIdx];
    return end - start;
}

int64_t SuffixTree::Node::edgeEnd(int64_t currentGlobalEnd, const std::vector<int64_t> &nodeEnds) const noexcept {
    return (endIdx == kLeafEnd) ? currentGlobalEnd : nodeEnds[endIdx];
}

// ─────────────────────────────────────────────────────────────────────────────
// textChar — returns the character at pos, or sentinel (255) at pos==_textLength
// ─────────────────────────────────────────────────────────────────────────────

static unsigned char textChar(const char *data, int64_t textLen, int64_t pos) noexcept {
    return (pos == textLen) ? static_cast<unsigned char>(255) : static_cast<unsigned char>(data[pos]);
}

// ─────────────────────────────────────────────────────────────────────────────
// Node allocation
// ─────────────────────────────────────────────────────────────────────────────

int64_t SuffixTree::newLeaf(int64_t start, int64_t suffixIndex) {
    // Leaf edges are open-ended and always point to the shared leaf sentinel.
    Node leaf;
    leaf.start = start;
    leaf.endIdx = kLeafEnd;

    // Leaf nodes represent a concrete suffix, so store its starting index now.
    leaf.suffixIndex = suffixIndex;

    // Newly created nodes start with no suffix link and no children.
    leaf.suffixLink = kRoot;
    leaf.children.fill(kInvalidIndex);

    _nodes.push_back(leaf);
    return static_cast<int64_t>(_nodes.size()) - 1;
}

int64_t SuffixTree::newInternal(int64_t start, int64_t end) {
    // Internal nodes own a dedicated end value stored in the end table.
    const int64_t endIdx = static_cast<int64_t>(_nodeEnds.size());
    _nodeEnds.push_back(end);

    Node node;
    node.start = start;
    node.endIdx = endIdx;

    // Internal nodes do not correspond to a suffix, so keep this invalid.
    node.suffixIndex = kInvalidIndex;

    // Internal nodes also begin with no suffix link and no children.
    node.suffixLink = kRoot;
    node.children.fill(kInvalidIndex);

    _nodes.push_back(node);
    return static_cast<int64_t>(_nodes.size()) - 1;
}

// ─────────────────────────────────────────────────────────────────────────────
// extendTree — one Ukkonen phase
// ─────────────────────────────────────────────────────────────────────────────

void SuffixTree::extendTree(int64_t pos) {
    _currentGlobalEnd = pos + 1;
    ++_remaining;

    int64_t lastNewInternal = kInvalidIndex;

    while (_remaining > 0) {
        // When the active length is zero, the active slot is determined by the current character.
        if (_activeLength == 0)
            _activeSlot = charToSlot(textChar(_textData, _textLength, pos));

        int64_t childIdx = _nodes[_activeNode].children[_activeSlot];

        if (childIdx == kInvalidIndex) {
            // Rule 2a: no edge exists, so create a new leaf for the current suffix.
            const int64_t suffixIdx = pos - _remaining + 1;
            const int64_t leaf = newLeaf(pos, suffixIdx);
            _nodes[_activeNode].children[_activeSlot] = leaf;

            // Resolve the suffix link chain for the last created internal node.
            if (lastNewInternal != kInvalidIndex) {
                _nodes[lastNewInternal].suffixLink = _activeNode;
                lastNewInternal = kInvalidIndex;
            }
        } else {
            // If the active point already covers the whole edge, move down to the child node.
            const int64_t edgeLen = _nodes[childIdx].edgeLength(_currentGlobalEnd, _nodeEnds);
            if (_activeLength >= edgeLen) {
                _activeSlot = charToSlot(textChar(_textData, _textLength,
                                                  pos - _remaining + 1 + edgeLen));
                _activeLength -= edgeLen;
                _activeNode = childIdx;
                continue;
            }

            // Compare the next character on the edge with the incoming character.
            const int64_t nextPos = _nodes[childIdx].start + _activeLength;
            const unsigned char nextOnEdge = textChar(_textData, _textLength, nextPos);
            const unsigned char incoming = textChar(_textData, _textLength, pos);

            if (nextOnEdge == incoming) {
                // Rule 3: the extension is already present, so increase the active length and stop.
                ++_activeLength;
                if (lastNewInternal != kInvalidIndex)
                    _nodes[lastNewInternal].suffixLink = _activeNode;
                break;
            }

            // split the edge and insert a new internal node plus a new leaf.
            const int64_t splitNode = newInternal(_nodes[childIdx].start, nextPos);
            _nodes[childIdx].start = nextPos;

            const int oldSlot = charToSlot(nextOnEdge);
            const int newSlot = charToSlot(incoming);
            _nodes[splitNode].children[oldSlot] = childIdx;

            const int64_t splitSuffixIdx = pos - _remaining + 1;
            const int64_t newLeafIdx = newLeaf(pos, splitSuffixIdx);
            _nodes[splitNode].children[newSlot] = newLeafIdx;

            _nodes[_activeNode].children[_activeSlot] = splitNode;

            // Link the previously created internal node to the newly split node.
            if (lastNewInternal != kInvalidIndex)
                _nodes[lastNewInternal].suffixLink = splitNode;
            lastNewInternal = splitNode;
        }

        --_remaining;

        // Recompute the active point for the next extension.
        if (_activeNode == kRoot && _activeLength > 0) {
            --_activeLength;
            const int64_t nextEdgePos = pos - _remaining + 1;
            _activeSlot = charToSlot(textChar(_textData, _textLength, nextEdgePos));
        } else {
            const int64_t link = _nodes[_activeNode].suffixLink;
            _activeNode = (link != kInvalidIndex) ? link : kRoot;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// initializeBuildState
// ─────────────────────────────────────────────────────────────────────────────

void SuffixTree::initializeBuildState() {
    // Start from a clean tree state before rebuilding.
    _nodes.clear();
    _nodeEnds.clear();

    // Reserve slot 0 for the root's dummy end value.
    _nodeEnds.push_back(0);

    // Initialize the root node with neutral/default suffix-tree values.
    Node rootNode;
    rootNode.start = 0;
    rootNode.endIdx = 0;
    rootNode.suffixIndex = kInvalidIndex;
    rootNode.suffixLink = kInvalidIndex;
    rootNode.children.fill(kInvalidIndex);
    _nodes.push_back(rootNode);

    // Reset the active point used by Ukkonen's algorithm.
    _activeNode = kRoot;
    _activeSlot = -1;
    _activeLength = 0;

    // No pending suffixes and no global end yet.
    _remaining = 0;
    _currentGlobalEnd = 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// buildSuffixTree
// ─────────────────────────────────────────────────────────────────────────────

void SuffixTree::buildSuffixTree() {
    if (_textLength == 0)
        return;

    initializeBuildState();

    const int64_t textLengthWithSentinel = _textLength + 1;

    // 2n+2 worst-case nodes: n leaves + (n-1) internal + root + sentinel leaf.
    _nodes.reserve(static_cast<size_t>(2 * textLengthWithSentinel + 4));
    _nodeEnds.reserve(static_cast<size_t>(textLengthWithSentinel + 2));

    for (int64_t pos = 0; pos < textLengthWithSentinel; ++pos)
        extendTree(pos);
}

// ─────────────────────────────────────────────────────────────────────────────
// collectLeaves — iterative DFS
// ─────────────────────────────────────────────────────────────────────────────

void SuffixTree::collectLeaves(int64_t subtreeRoot,
                               std::vector<STSearchResult> &out,
                               int64_t patternLength) const {
    std::vector<int64_t> stack;
    stack.reserve(64);
    stack.push_back(subtreeRoot);

    while (!stack.empty()) {
        const int64_t idx = stack.back();
        stack.pop_back();

        const Node &node = _nodes[idx];

        if (node.suffixIndex != kInvalidIndex) {
            // Leaf: suffixIndex holds the 0-based start position.
            if (node.suffixIndex < _textLength)
                out.push_back({node.suffixIndex, patternLength});
            continue; // leaves have no children
        }

        // Internal node: enqueue children (reverse order for left-to-right DFS).
        for (int slot = kAlphabetSize - 1; slot >= 0; --slot) {
            const int64_t child = node.children[slot];
            if (child != kInvalidIndex)
                stack.push_back(child);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// search — O(m) descent + DFS leaf collection
// ─────────────────────────────────────────────────────────────────────────────

std::vector<STSearchResult>
SuffixTree::search(const std::string &pattern) const {
    // If the tree is not built yet or the pattern is empty, there can be no matches.
    if (!_isReady || pattern.empty())
        return {};

    const int64_t patternLength = static_cast<int64_t>(pattern.size());
    int64_t currentNodeIndex = kRoot;
    int64_t patternPos = 0;

    // Collect all matches in the subtree below the current node,
    // then sort them by offset for stable, predictable output.
    const auto collectSortedResults = [&](int64_t subtreeRoot) {
        std::vector<STSearchResult> results;
        collectLeaves(subtreeRoot, results, patternLength);
        std::sort(results.begin(), results.end(),
                  [](const STSearchResult &a, const STSearchResult &b) {
                      return a.offset < b.offset;
                  });
        return results;
    };

    // Traverse the tree while characters from the pattern still need to be matched.
    while (patternPos < patternLength) {
        const unsigned char patternChar = static_cast<unsigned char>(pattern[patternPos]);
        const int slot = charToSlot(patternChar);
        const int64_t childIndex = _nodes[currentNodeIndex].children[slot];

        // No outgoing edge for this character means the pattern is absent.
        if (childIndex == kInvalidIndex)
            return {};

        const Node &child = _nodes[childIndex];
        const int64_t edgeStart = child.start;
        const int64_t edgeEnd = child.edgeEnd(_currentGlobalEnd, _nodeEnds);

        // Compare the pattern against the current edge label character by character.
        for (int64_t edgePos = edgeStart;
             edgePos < edgeEnd && patternPos < patternLength;
             ++edgePos, ++patternPos) {
            if (textChar(_textData, _textLength, edgePos) !=
                static_cast<unsigned char>(pattern[patternPos]))
                return {};
        }

        // If the pattern is not fully matched yet, continue descending into the child node.
        if (patternPos < patternLength) {
            currentNodeIndex = childIndex;
            continue;
        }

        // Pattern fully matched inside this subtree; collect all leaf matches below it.
        return collectSortedResults(childIndex);
    }

    // If the loop ends exactly at a node boundary, collect matches from that node's subtree.
    return collectSortedResults(currentNodeIndex);
}