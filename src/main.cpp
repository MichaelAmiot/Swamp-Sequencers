#include "SwampSeqLib/genome_mapper.h"
#include "SwampSeqLib/suffix_array.h"
#include <iostream>
#include <string_view>

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

  mappedData[1000000] = '\0';

  SuffixArray SA(mappedData);
  std::cout << "Suffix array built, searching...\n";
  auto res = SA.search("ACT");
  std::cout << res.size() << " Patterns found at:\n";
  for (auto pair : res) {
    std::cout << "Offset: " << pair.offset << "\n";
  }
}
