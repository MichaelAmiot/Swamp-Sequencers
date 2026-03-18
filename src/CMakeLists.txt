# Create core functionality library for testing
add_library(SwampSeqLib STATIC "")
target_sources(SwampSeqLib PRIVATE
  suffix_array.cpp
  suffix_tree.cpp
)

target_include_directories(SwampSeqLib PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_compile_options(SwampSeqLib PRIVATE -Wall -Wextra -Wpedantic)

# Create main executable
add_executable(SwampSequencer "")
target_sources(SwampSequencer PRIVATE main.cpp)

# Link the main executable to the library
target_link_libraries(SwampSequencer PRIVATE SwampSeqLib)



