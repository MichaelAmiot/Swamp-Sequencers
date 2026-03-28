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

# Genome Download
FetchContent_Declare(
    ecoli_genome
    URL "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000008865.2/download?include_annotation_type=GENOME_FASTA"
    DOWNLOAD_NO_EXTRACT FALSE  # NCBI sends a ZIP, set to TRUE if you want the raw zip
)

FetchContent_MakeAvailable(ecoli_genome)

# Locate the downloaded file
set(GENOME_SOURCE_FILE "${ecoli_genome_SOURCE_DIR}/ncbi_dataset/data/GCF_000008865.2/GCF_000008865.2_ASM886v2_genomic.fna")

# Copy to Binary Directory
# This ensures that 'ecoli.fna' is always sitting next to SwampSequencer.exe
add_custom_command(
    TARGET SwampSequencer POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${GENOME_SOURCE_FILE}
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ecoli.fna
    COMMENT "Copying E. coli genome to binary directory..."
)

# Link the main executable to the library
target_link_libraries(SwampSequencer PRIVATE SwampSeqLib)
target_link_libraries(SwampSequencer PRIVATE 
  ftxui::screen
  ftxui::dom
  ftxui::component
)



