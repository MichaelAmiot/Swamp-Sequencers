#include "SwampSeqLib/genome_mapper.h"
#include "SwampSeqLib/suffix_array.h"
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <path>\n";
    return 1;
  }

  GenomeMapper data(argv[1]);
  if (!data.isValid()) {
    std::cerr << "Could not open: " << argv[1] << "\n";
    return 1;
  }
  char *mappedData = data.data();
  std::cout << "File size: " << data.size() << "\n" << std::flush;
  std::cout << "Building suffix array.\n";
  SuffixArray SA(data);
  std::cout << "Suffix array built, searching...\n";
  auto res =
      SA.search("CATGATTAGGAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
  std::cout << "pattern found at:\n";
  for (auto pair : res) {
    std::cout << "Offset: " << pair.offset << "\n";
  }
}
