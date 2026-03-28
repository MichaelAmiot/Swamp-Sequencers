# Set minimum CMake version required
cmake_minimum_required(VERSION 3.15)

# Define project
project(SwampSequencer VERSION 1.0 LANGUAGES CXX)

# Define C++ Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Centralize binary output for easier testing
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# Enable CTest
enable_testing()

# Fetch GTest framework
include(fetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.17.0
)
FetchContent_Declare(ftxui
  GIT_REPOSITORY https://github.com/ArthurSonzogni/ftxui
  GIT_TAG v6.1.9
)
FetchContent_MakeAvailable(ftxui)




# Add src subdirectory
add_subdirectory(src)

# Option for Windows compatibility
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_subdirectory(tests)
