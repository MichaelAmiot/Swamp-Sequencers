#include "SwampSeqLib/genome_mapper.h"
#include <cstddef>
#include <fileapi.h>
#include <handleapi.h>
#include <memoryapi.h>
#include <string>
#include <winnt.h>

GenomeMapper::GenomeMapper(const std::string &filePath) {
#ifdef OS_WINDOWS
  // Open file
  _fileHandle = CreateFileA(filePath.c_str(), GENERIC_READ, FILE_SHARE_READ,
                            NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
  if (_fileHandle == INVALID_HANDLE_VALUE)
    return;

  // Get size
  LARGE_INTEGER size;
  if (!GetFileSizeEx(_fileHandle, &size))
    return;
  _fileSize = size.QuadPart;

  // Get map handle
  _memoryHandle =
      CreateFileMappingA(_fileHandle, NULL, PAGE_READONLY, 0, 0, NULL);
  if (_memoryHandle == NULL)
    return;

  // Get mapped data
  _mappedData = MapViewOfFile(_memoryHandle, FILE_MAP_READ, 0, 0, 0);
  if (_mappedData != NULL)
    _isValid = true;

#elif defined(OS_POSIX)
  // Open file
  _fileDescriptor = open(filePath.c_str(), O_RDONLY);
  if (_fileDescriptor == -1)
    return;

  // Get size
  struct stat sb;
  if (fstat(_fileDescriptor, &sb) == -1)
    return;
  _fileSize = sb.st_size;

  // Map data
  _mappedData =
      mmap(nullptr, _fileSize, PROT_READ, MAP_PRIVATE, _fileDescriptor, 0);
  if (_mappedData != MAP_FAILED)
    _isValid = true;
#endif
}

GenomeMapper::~GenomeMapper() {
  // Handle clean-up based on OS
#ifdef OS_WINDOWS
  if (_mappedData)
    UnmapViewOfFile(_mappedData);
  if (_memoryHandle)
    CloseHandle(_memoryHandle);
  if (_fileHandle != INVALID_HANDLE_VALUE)
    CloseHandle(_fileHandle);
#elif defined(OS_POSIX)
  if (_isValid)
    munmap(_mappedData, _fileSize);
  if (_fileDescriptor != -1)
    close(_fileDescriptor);
#endif
}

const char *GenomeMapper::data() const {
  return static_cast<const char *>(_mappedData);
}
size_t GenomeMapper::size() const { return _fileSize; }
bool GenomeMapper::isValid() const { return _isValid; }
