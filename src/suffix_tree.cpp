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
// kAlphabetSize
// ─────────────────────────────────────────────────────────────────────────────

// Map a raw byte to a children[] slot index.
// Sentinel (-1 cast to unsigned = 255) maps to slot 5.
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

    constexpr size_t MAX_BP = 500'000'000;
    if (mapper.size() > MAX_BP)
        throw std::runtime_error(
            "SuffixTree: genome exceeds 500 Mbp hard cap. "
            "Use SuffixArray for large genomes.");

    _textData = mapper.data();
    _textLength = static_cast<uint32_t>(mapper.size());

    buildSuffixTree();
    _isReady = true;
}

SuffixTree::SuffixTree(const std::string &text) {
    constexpr size_t MAX_BP = 500'000'000;
    if (text.size() > MAX_BP)
        throw std::runtime_error(
            "SuffixTree: text exceeds 500 Mbp hard cap. "
            "Use SuffixArray for large genomes.");

    _owned = text;
    _textData = _owned.c_str();
    _textLength = static_cast<uint32_t>(_owned.size());

    buildSuffixTree();
    _isReady = true;
}

// ─────────────────────────────────────────────────────────────────────────────
// Ukkonen's Algorithm — O(n) time, O(n) space
//
// Node memory: 40 bytes flat (6 × uint32_t children + 4 × uint32_t fields).
// No heap allocation per node — all nodes live in _nodes vector.
// Internal node end values live in _nodeEnds; leaves share _currentGlobalEnd via
// the kLeafEndvsentinel in endIdx.
// ─────────────────────────────────────────────────────────────────────────────

uint32_t SuffixTree::newLeaf(uint32_t start) {
    Node leaf;
    leaf.start = start;
    leaf.endIdx = kLeafEnd; // open end — reads _currentGlobalEnd at query time
    leaf.suffixIndex = kInvalidIndex; // filled by annotateSuffixIndices()
    leaf.suffixLink = kRoot;
    _nodes.push_back(leaf);
    return static_cast<uint32_t>(_nodes.size()) - 1;
}

uint32_t SuffixTree::newInternal(uint32_t start, uint32_t end) {
    const uint32_t endIdx = static_cast<uint32_t>(_nodeEnds.size());
    _nodeEnds.push_back(end);

    Node node;
    node.start = start;
    node.endIdx = endIdx;
    node.suffixIndex = kInvalidIndex;
    node.suffixLink = kRoot;
    _nodes.push_back(node);
    return static_cast<uint32_t>(_nodes.size()) - 1;
}

// ── Helper: character at position pos in the text (sentinel at pos == _textLength) ──

// Edge length — for leaves *end is _currentGlobalEnd at query time.
uint32_t SuffixTree::Node::edgeLength(uint32_t currentGlobalEnd, const std::vector<uint32_t> &nodeEnds) const noexcept {
        const uint32_t edgeEnd = (endIdx == kLeafEnd) ? currentGlobalEnd : nodeEnds[endIdx];
        return edgeEnd - start;
}

uint32_t SuffixTree::Node::edgeEnd(uint32_t globalEnd, const std::vector<uint32_t> &nodeEnds) const noexcept {
    const uint32_t resolvedEnd = (endIdx == kLeafEnd) ? globalEnd : nodeEnds[endIdx];
    return resolvedEnd;
}

// Forward-declare the free helper used inside extendTree.
static unsigned char textChar(const char *data, uint32_t num, uint32_t pos) {
    return (pos == num) ? static_cast<unsigned char>(255) : static_cast<unsigned char>(data[pos]); // sentinel byte
}

// ─────────────────────────────────────────────────────────────────────────────
// extendTree — one Ukkonen phase
// ─────────────────────────────────────────────────────────────────────────────

void SuffixTree::extendTree(uint32_t pos) {
    // Extend the implicit suffix tree with the next character.
    _currentGlobalEnd = pos + 1;
    ++_remaining;

    // Tracks the most recently created internal node so we can wire suffix links.
    uint32_t lastNewInternal = kInvalidIndex;

    // Process all pending suffix extensions for this phase.
    while (_remaining > 0) {
        // If we are at a node, the active edge starts with the current input character.
        if (_activeLength == 0) {
            _activeSlot = SuffixTree::charToSlot(textChar(_textData, _textLength, pos));
        }

        // Child reached from the active node via the active edge.
        uint32_t childIdx = _nodes[_activeNode].children[_activeSlot];

        if (childIdx == kInvalidIndex) {
            // Rule 2a: no matching edge exists, so create a new leaf here.
            const uint32_t leaf = newLeaf(pos);
            _nodes[_activeNode].children[_activeSlot] = leaf;

            // If we created an internal node in the previous extension,
            // link it to the current active node.
            if (lastNewInternal != kInvalidIndex) {
                _nodes[lastNewInternal].suffixLink = _activeNode;
                lastNewInternal = kInvalidIndex;
            }
        } else {
            // If the active point spans the full edge, move down to the child node.
            const uint32_t edgeLen = _nodes[childIdx].edgeLength(_currentGlobalEnd, _nodeEnds);
            if (_activeLength >= edgeLen) {
                _activeSlot = SuffixTree::charToSlot(textChar(_textData, _textLength, pos - _remaining + 1 + edgeLen));
                _activeLength -= edgeLen;
                _activeNode = childIdx;
                continue;
            }

            // Compare the next character on the edge with the incoming character.
            const uint32_t nextPos = _nodes[childIdx].start + _activeLength;
            const unsigned char nextOnEdge = textChar(_textData, _textLength, nextPos);
            const unsigned char incoming = textChar(_textData, _textLength, pos);

            if (nextOnEdge == incoming) {
                // Rule 3: the extension is already present, so this phase stops.
                ++_activeLength;

                // Any pending suffix-link chain ends at the current active node.
                if (lastNewInternal != kInvalidIndex)
                    _nodes[lastNewInternal].suffixLink = _activeNode;

                break;
            }

            // Rule 2b: mismatch found, so split the edge and create a new branch.
            const uint32_t splitNode = newInternal(_nodes[childIdx].start, nextPos);

            // The old child edge now starts after the split point.
            _nodes[childIdx].start = nextPos;

            // Attach both the old child and the new leaf to the split node.
            const int oldSlot = charToSlot(nextOnEdge);
            const int newSlot = charToSlot(incoming);
            _nodes[splitNode].children[oldSlot] = childIdx;

            const uint32_t newLeafIdx = newLeaf(pos);
            _nodes[splitNode].children[newSlot] = newLeafIdx;

            // Replace the original child with the split node.
            _nodes[_activeNode].children[_activeSlot] = splitNode;

            // Maintain the suffix-link chain between newly created internal nodes.
            if (lastNewInternal != kInvalidIndex)
                _nodes[lastNewInternal].suffixLink = splitNode;
            lastNewInternal = splitNode;
        }

        // One pending suffix extension has been resolved.
        --_remaining;

        // Update the active point for the next extension.
        if (_activeNode == kRoot
            && _activeLength > 0) {
            --_activeLength;
            const uint32_t nextEdgePos = pos - _remaining + 1;
            _activeSlot = SuffixTree::charToSlot(textChar(_textData, _textLength, nextEdgePos));
        } else {
            const uint32_t link = _nodes[_activeNode].suffixLink;
            _activeNode = (link != kInvalidIndex) ? link : kRoot;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// annotateSuffixIndices — iterative DFS
//
// suffixIndex for leaf = (_textLength+ 1) - labelHeight
// (_textLength+ 1 because the sentinel adds one to every kRoot-to-leaf path length.)
// ─────────────────────────────────────────────────────────────────────────────

void SuffixTree::annotateSuffixIndices() {
    // Stack frame for iterative DFS:
    // nodeIdx - current node being processed
    // height  - total path length from kRoot to this node
    struct Frame {
        uint32_t nodeIdx;
        uint32_t height;
    };

    // Start traversal from the kRoot with zero path height.
    std::vector<Frame> stack;
    stack.reserve(1024);
    stack.push_back({kRoot, 0u});

    while (!stack.empty()) {
        // Process the most recently discovered node.
        auto [nodeIdx, height] = stack.back();
        stack.pop_back();

        Node &node = _nodes[nodeIdx];
        bool isLeaf = true;

        // Visit all children and accumulate path length along each edge.
        for (int slot = 0; slot < kAlphabetSize; ++slot) {
            const uint32_t child = node.children[slot];
            if (child == kInvalidIndex) {
                continue;
            }

            isLeaf = false;

            // Child height = current path height + edge length to the child.
            const uint32_t childHeight = height + _nodes[child].edgeLength(_currentGlobalEnd, _nodeEnds);
            stack.push_back({child, childHeight});
        }

        // A leaf suffix corresponds to the text position:
        //   suffixIndex = (text length including sentinel) - path height.
        if (isLeaf) {
            node.suffixIndex = (_textLength + 1) - height;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// buildSuffixTree 
// ─────────────────────────────────────────────────────────────────────────────

void SuffixTree::initializeBuildState() {
    // Reset all mutable build state before starting a fresh construction.
    _nodes.clear();
    _nodeEnds.clear();

    // Reserve a private end value for the root node to keep indexing consistent.
    _nodeEnds.push_back(0);

    // Create the root node explicitly and store it at index kRoot (0).
    Node rootNode;
    rootNode.start = 0;
    rootNode.endIdx = 0;
    rootNode.suffixIndex = kInvalidIndex;
    rootNode.suffixLink = kInvalidIndex;
    _nodes.push_back(rootNode);

    // Ukkonen's active point starts at the root with no pending extension.
    _activeNode = kRoot;
    _activeSlot = -1;
    _activeLength = 0;
    _remaining = 0;
    _currentGlobalEnd = 0;
}

void SuffixTree::buildSuffixTree() {
    if (_textLength == 0) {
        return;
    }

    // Build the tree for the text plus one sentinel character.
    const uint32_t textLengthWithSentinel = _textLength + 1;

    // Reserve enough space for the worst-case number of nodes and end values.
    const uint32_t maxNodeCount = 2 * textLengthWithSentinel + 4;
    _nodes.reserve(maxNodeCount);
    _nodeEnds.reserve(textLengthWithSentinel + 2);

    // Initialize the root node and all build-time variables.
    initializeBuildState();

    // Extend the tree one character at a time, including the sentinel.
    for (uint32_t pos = 0; pos < textLengthWithSentinel; ++pos) {
        extendTree(pos);
    }

    // Assign suffix indices to leaves after the structure is fully built.
    annotateSuffixIndices();
}

// ─────────────────────────────────────────────────────────────────────────────
// collectLeaves — iterative DFS
// ─────────────────────────────────────────────────────────────────────────────

void SuffixTree::collectLeaves(uint32_t subtreeRoot, std::vector<STSearchResult> &out, uint32_t patternLength) const {
    // Use an explicit stack instead of recursion to avoid deep call chains.
    static constexpr std::size_t INITIAL_STACK_CAPACITY = 1024;

    std::vector<uint32_t> stack;
    stack.reserve(INITIAL_STACK_CAPACITY);
    stack.push_back(subtreeRoot);

    while (!stack.empty()) {
        // Process the next node in depth-first order.
        const uint32_t currentNodeIndex = stack.back();
        stack.pop_back();

        const Node &currentNode = _nodes[currentNodeIndex];
        bool hasChildren = false;

        // Push all existing children onto the stack.
        for (std::size_t childSlot = 0;
             childSlot < static_cast<std::size_t>(kAlphabetSize); ++childSlot) {
            const uint32_t childIndex = currentNode.children[childSlot];
            if (childIndex == kInvalidIndex) {
                continue;
            }

            hasChildren = true;
            stack.push_back(childIndex);
        }

        // A valid result is a leaf with a suffix index inside the text bounds.
        const bool isValidLeaf = !hasChildren && currentNode.suffixIndex != kInvalidIndex && currentNode.suffixIndex <
                                 _textLength;

        if (isValidLeaf) {
            out.push_back({currentNode.suffixIndex, patternLength});
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// search — O(m) descent + DFS leaf collection
// ─────────────────────────────────────────────────────────────────────────────

std::vector<STSearchResult>
SuffixTree::search(const std::string &pattern) const {
    // Fast fail: tree is not built yet or the query is empty.
    if (!_isReady || pattern.empty()) {
        return {};
    }

    const uint32_t patternLength = static_cast<uint32_t>(pattern.size());
    uint32_t currentNodeIndex = kRoot;
    uint32_t patternPos = 0;

    // Collect all leaves under subtreeRoot and sort results by offset.
    const auto collectSortedResults = [&](uint32_t subtreeRoot) {
        std::vector<STSearchResult> results;
        collectLeaves(subtreeRoot, results, patternLength);

        std::sort(results.begin(), results.end(),
                  [](const STSearchResult &lhs, const STSearchResult &rhs) {
                      return lhs.offset < rhs.offset;});

        return results;
    };

    // Walk down the suffix tree while matching the pattern against edge labels.
    while (patternPos < patternLength) {
        const unsigned char patternChar = static_cast<unsigned char>(pattern[patternPos]);
        const int slot = charToSlot(patternChar);
        const uint32_t childIndex = _nodes[currentNodeIndex].children[slot];

        // No outgoing edge for the next pattern character means no match.
        if (childIndex == kInvalidIndex) {
            return {};
        }

        const Node &child = _nodes[childIndex];
        const uint32_t edgeStart = child.start;
        const uint32_t edgeEnd = child.edgeEnd(_currentGlobalEnd, _nodeEnds);

        // Compare pattern characters with the current edge label.
        for (uint32_t edgePos = edgeStart;
             edgePos < edgeEnd && patternPos < patternLength;
             ++edgePos, ++patternPos) {
            if (textChar(_textData, _textLength, edgePos) !=
                static_cast<unsigned char>(pattern[patternPos])) {
                return {};
            }
        }

        // Pattern still has characters left: continue from the child node.
        if (patternPos < patternLength) {
            currentNodeIndex = childIndex;
            continue;
        }

        // Pattern ended inside this subtree; collect all matching leaves.
        return collectSortedResults(childIndex);
    }

    // Pattern ended exactly at an internal node.
    return collectSortedResults(currentNodeIndex);
}