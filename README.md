# SwampSequencer
## Overview
Given a massive DNA genome database, finding where a specific gene sequence pattern appears is the core operation behind disease detection, ancestry matching, and drug research.

This project benchmarks the performance, memory efficiency, and build constraints of two O(n) string-matching data structures: the Suffix Array (using the SA-IS algorithm) and the Suffix Tree (using Ukkonen's algorithm). We test these structures against real human genome data, specifically NCBI GenBank's Human Chromosome 22 (GRCh38), which contains over 100 million base pairs.
## Getting Started
### 1. Clone the repository
```Bash
git clone https://github.com/MichaelAmiot/Swamp-Sequencers
cd Swamp-Sequencers
```

### 2. Obtain the data
#### Download the data
The project uses a large dataset which doesn't fit in a GitHub repo, but can be obtained [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) by clicking "Download", selecting "Genome Sequences (FASTA)", and clicking Download again.
#### Moving the data
After downloading the data, create a data directory in the Swamp-Sequencers project directory.
```Bash
# Swamp-Sequencers/
mkdir data
```
Then, place the downloaded zip file ("ncbi_dataset.zip" unless you've changed it) into the data directory and unzip the file.
### 3. Build the project
The project uses an out-of-source build.
```Bash
mkdir build
cd build
cmake ..
make
```
### Usage
Run the CLI interface to query the dataset
```Bash
./bin/SwampSequencer
```
Or run the automated Google Test suite
```Bash
./bin/SwampTests
```
## Attribution
    - Michael Amiot
    - Jack Hendrix
    - Sebastian Mejia
