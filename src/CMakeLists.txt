cmake_minimum_required(VERSION 3.20)
project(SwampSequencers CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create core functionality library for testing
add_library(SwampSeqLib STATIC "")
target_sources(SwampSeqLib PRIVATE
  suffix_array.cpp
  suffix_tree.cpp
  genome_mapper.cpp
)

target_include_directories(SwampSeqLib PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_compile_options(SwampSeqLib PRIVATE -Wall -Wextra -Wpedantic)

# Create main executable
add_executable(SwampSequencer "")
target_sources(SwampSequencer PRIVATE main.cpp)

# Link the main executable to the library
target_link_libraries(SwampSequencer PRIVATE SwampSeqLib)



