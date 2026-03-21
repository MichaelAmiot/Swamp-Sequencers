#include "SwampSeqLib/genome_mapper.h"
#include <iostream>
int main(int argc, char *argv[]) {
  GenomeMapper mapped(argv[1]);
  const char *data = mapped.data();
  if (!mapped.isValid()) {
    std::cerr << "Could not open file: " << argv[1] << "\n";
  }
  std::cout << "Size of data: " << mapped.size() << "\n";
  std::cout << "First 1000 characters:\n";
  for (int i = 0; i < 1000; i++) {
    std::cout << data[i];
  }
  return 0;
}
