#pragma once
#include <string>

// Define memory mapper based on OS
#if defined(_WIN32) || defined(_WIN64)
#define OS_WINDOWS
#define NOMINMAX
#include <windows.h>
#else
#define OS_POSIX
#include <fcntl.h> open()
#include <sys/mmap.h> mmap(), mumap()
#include <sys/stat.h> fstat()
#include <unistd.h>
#endif

class GenomeMapper {
#ifdef OS_WINDOWS
  HANDLE _fileHandle = INVALID_HANDLE_VALUE;
  HANDLE _memoryHandle = NULL;
#elif defined(OS_POSIX)
  int _fileDescriptor = -1;
#endif
  size_t _fileSize = 0;
  void *_mappedData = nullptr;
  bool _isValid = false;

public:
  // Constructor/Destructor
  GenomeMapper(const std::string &filepath);
  ~GenomeMapper();

  // Remove copying
  GenomeMapper(const GenomeMapper &) = delete;
  GenomeMapper operator=(const GenomeMapper &) = delete;

  // Accessors
  const char *data() const;
  size_t size() const;
  bool isValid() const;
};
