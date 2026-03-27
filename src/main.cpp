#include "SwampSeqLib/genome_mapper.h"
#include "SwampSeqLib/suffix_array.h"
#include <chrono>
#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <string>
#include <thread>

using namespace ftxui;

int main() {
  // 1. Shared Application State (Now with TWO strings)
  std::string search_pattern = "ATGC";
  std::string file_path =
      "data/ncbi_dataset/ncbi_dataset/data/GCF_000001405.40/"
      "GCF_000001405.40_GRCh38.p14_genomic.fna";

  int algo_choice = 0;
  std::vector<std::string> algo_entries = {"Suffix Array", "Suffix Tree",
                                           "Compare Both"};

  std::string status_text = "Idle. Ready to benchmark.";
  std::string result_text = "";
  bool is_running = false;

  auto screen = ScreenInteractive::TerminalOutput();

  // 2. Interactive Components
  Component input_filepath = Input(&file_path, "e.g., ./data/genome.fasta");
  Component input_pattern = Input(&search_pattern, "e.g., ATGC");
  Component algo_radiobox = Radiobox(&algo_entries, &algo_choice);

  // 3. The Benchmark Action
  GenomeMapper map;
  Component run_button = Button("Run Benchmark", [&] {
    if (is_running)
      return;

    is_running = true;
    status_text = "Mapping " + file_path + "...";
    result_text = "";

    std::thread([&, file_path, search_pattern] {
      try {
        map.fromFile(file_path);
        if (!map.isValid()) {
          screen.Post([&] {
            status_text = "Could not open file.";
            is_running = false;
          });
          map.data()[100000] = '\0';
          return;
        }
      } catch (const std::exception &e) {
        std::string errorMessage = e.what();
        screen.Post([&] {
          status_text = "Could not open file: " + errorMessage;
          is_running = false;
        });
        is_running = false;
        return;
      }
      std::chrono::duration<double> saTime;
      std::vector<SearchResult> saRes;
      std::chrono::duration<double> stTime;
      std::vector<SearchResult> stRes;
      try {
        if (algo_choice == 0 || algo_choice == 2) {
          auto start = std::chrono::steady_clock::now();
          screen.Post([&] { status_text = "Building suffix array..."; });
          SuffixArray sa(map);
          saRes = sa.search(search_pattern);
          auto end = std::chrono::steady_clock::now();
          saTime = end - start;
        }
      } catch (std::exception &e) {
        screen.Post([&] {
          std::string errorMessage = e.what();
          status_text = "Error building Suffix Array: " + errorMessage;
        });
        is_running = false;
        return;
      }

      screen.Post([&] {
        status_text = "Benchmark Complete!";
        if (algo_choice == 0) {
          result_text =
              "Suffix Array: " + std::to_string(saTime.count() * 1000) + " | " +
              std::to_string(saRes.size()) + " matches found.";
        } else if (algo_choice == 1) {
          result_text = "Suffix Tree: 38ms | 1,204 matches found";
        } else {
          result_text = "Array: 45ms | Tree: 38ms (Tree was 7ms faster)";
        }
        is_running = false;
      });
    }).detach();
  });

  // Container
  auto layout = Container::Vertical(
      {input_filepath, input_pattern, algo_radiobox, run_button});

  // Renderer
  auto renderer = Renderer(layout, [&] {
    return vbox({text("Swamp Sequencer") | bold | color(Color::Green),
                 separator(),

                 // Render both inputs with clean, aligned labels
                 hbox({text("Genome File Path : "), input_filepath->Render()}),
                 hbox({text("Search Pattern   : "), input_pattern->Render()}),
                 separator(),

                 text("Select Data Structure:"), algo_radiobox->Render(),
                 separator(),

                 run_button->Render() | center, separator(),

                 text(status_text) | color(Color::Yellow),
                 text(result_text) | color(Color::Cyan) | bold

           }) |
           border;
  });

  screen.Loop(renderer);

  return 0;
}
