# SwampSequencer
## Overview
Given a massive DNA genome database, finding where a specific gene sequence pattern appears is the core operation behind disease detection, ancestry matching, and drug research.

This project benchmarks the performance, memory efficiency, and build constraints of two O(n) string-matching data structures: the Suffix Array (using the SA-IS algorithm) and the Suffix Tree (using Ukkonen's algorithm). We test these structures against real genome data, specifically NCBI GenBank's E. Coli strain Sakai sample chromosome found [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000008865.2/).
## Getting Started
### 1. Clone the repository
```Bash
git clone https://github.com/MichaelAmiot/Swamp-Sequencers
cd Swamp-Sequencers
```

### 2. Build the project
The project uses an out-of-source build. The first time the project is built it will also download the E. Coli reference genome and place it next to the binary.
```Bash
# Swamp-Sequencers/
mkdir build
cd build
# Swamp-Sequencers/build
cmake ..
make
```
### Usage
Run the TUI interface to query the dataset.
```Bash
# Swamp-Sequencers/build
bin/SwampSequencer.exe
```
Or run the automated Google Test suite
```Bash
# Swamp-Sequencers/build
bin/SwampTests.exe
```
## Attribution
    - Michael Amiot
    - Jack Hendrix
    - Sebastian Mejia
