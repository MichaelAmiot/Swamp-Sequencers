#include "SwampSeqLib/genome_mapper.h"
#include "SwampSeqLib/suffix_array.h"
#include "SwampSeqLib/suffix_tree.h"
#include <algorithm>
#include <chrono>
#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/dom/elements.hpp>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace ftxui;

int main() {
  // Shared Application State
  std::string searchPattern = "ATGC";
  std::string filePath = "data/ecoli/ecoli.fna";
  int algoChoice = 0;
  std::vector<std::string> algoEntries = {"Suffix Array", "Suffix Tree",
                                          "Compare Both"};

  std::string statusText = "Idle. Ready to benchmark.";
  std::string resultText = "";
  std::vector<std::string> previewLines;
  int selectedResult = 0;
  bool isRunning = false;

  auto screen = ScreenInteractive::TerminalOutput();

  // Interactive Components
  Component inputFilepath = Input(&filePath, "e.g., ./data/genome.fna");
  Component inputPattern = Input(&searchPattern, "e.g., ATGC");
  Component algoRadiobox = Radiobox(&algoEntries, &algoChoice);

  // Menu component for high-performance scrolling
  Component resultList = Menu(&previewLines, &selectedResult);

  // Benchmarking Logic
  Component runButton = Button("Run Benchmark", [&] {
    if (isRunning)
      return;

    isRunning = true;
    statusText = "Mapping " + filePath + "...";
    resultText = "";
    previewLines.clear();

    std::thread([&, filePath, searchPattern] {
      GenomeMapper map;
      try {
        map.fromFile(filePath);
        if (!map.isValid()) {
          screen.Post([&] {
            statusText = "Could not open file.";
            isRunning = false;
          });
          screen.Post(Event::Custom);
          return;
        }
      } catch (const std::exception &e) {
        std::string errMsg = e.what();
        screen.Post([&, errMsg] {
          statusText = "Error: " + errMsg;
          isRunning = false;
        });
        screen.Post(Event::Custom);
        return;
      }

      // Suffix Array Search
      std::vector<SearchResult> saRes;
      double saTimeMs;
      int saMatchCount;
      if (algoChoice == 0 || algoChoice == 2) {
        SuffixArray sa(map);

        auto saStart = std::chrono::steady_clock::now();
        saRes = sa.search(searchPattern);
        auto saEnd = std::chrono::steady_clock::now();

        saMatchCount = static_cast<int>(saRes.size());
        saTimeMs =
            std::chrono::duration<double>(saEnd - saStart).count() * 1000.0;
      }
      // Suffix Tree Search
      std::vector<STSearchResult> stRes;
      double stTimeMs;
      int stMatchCount;
      if (algoChoice == 1 || algoChoice == 2) {
        SuffixTree st(map);

        auto stStart = std::chrono::steady_clock::now();
        stRes = st.search(searchPattern);
        auto stEnd = std::chrono::steady_clock::now();

        stMatchCount = static_cast<int>(saRes.size());
        stTimeMs =
            std::chrono::duration<double>(stEnd - stStart).count() * 1000.0;
      }

      // Build previews on background thread
      std::vector<std::string> localStrings;
      const char *data = map.data();
      int previewLimit = std::min<int>(saMatchCount, 2000);
      int contextLen = 25;

      for (int i = 0; i < previewLimit; i++) {
        size_t idx = saRes[i].offset;
        size_t startIdx = (idx > contextLen) ? (idx - contextLen) : 0;
        size_t endIdx = std::min<size_t>(
            map.size(), idx + searchPattern.size() + contextLen);

        std::string snippet(data + startIdx, endIdx - startIdx);
        localStrings.push_back("Offset(" + std::to_string(idx) + "): ..." +
                               snippet + "...");
      }

      // Format the time to 1 decimal place
      std::stringstream ss;
      ss.imbue(std::locale::classic());
      ss << std::fixed << std::setprecision(1) << saTimeMs;
      std::string saFormattedTime = ss.str();

      ss.str("");
      ss.clear();
      ss << std::fixed << std::setprecision(1) << stTimeMs;
      std::string stFormattedTime = ss.str();

      // Update UI state
      screen.Post([&, saMatchCount, saFormattedTime, stMatchCount,
                   stFormattedTime, localStrings] {
        previewLines = std::move(localStrings);
        statusText = "Benchmark Complete!";
        resultText = "";
        if (algoChoice == 0 || algoChoice == 2) {
          resultText += "Suffix Array: " + saFormattedTime + "ms, " +
                        std::to_string(saMatchCount) + " matches found.";
        }
        if (algoChoice == 2) {
          resultText += " | ";
        }
        if (algoChoice == 1 || algoChoice == 2) {
          resultText += "Suffix Tree: " + stFormattedTime + "ms, " +
                        std::to_string(stMatchCount) + " matches found.";
        }
        isRunning = false;
      });

      // Trigger redraw immediately
      screen.Post(Event::Custom);
    }).detach();
  });

  // Layout & Navigation Container
  auto layout = Container::Vertical({
      inputFilepath,
      inputPattern,
      algoRadiobox,
      runButton,
      resultList,
  });

  // Renderer
  auto renderer = Renderer(layout, [&] {
    return vbox({text("SWAMP SEQUENCER") | bold | color(Color::Green) | center,
                 separator(),

                 hbox({text("Genome File Path : "), inputFilepath->Render()}),
                 hbox({text("Search Pattern   : "), inputPattern->Render()}),
                 separator(),

                 text("Select Data Structure:"), algoRadiobox->Render(),
                 separator(),

                 runButton->Render() | center, separator(),

                 text(statusText) | color(Color::Yellow),
                 text(resultText) | color(Color::Cyan) | bold,

                 // Results window
                 vbox({text("Results (Arrows to scroll)") | dim, separator(),
                       resultList->Render() | vscroll_indicator | frame |
                           size(HEIGHT, LESS_THAN, 15)}) |
                     border | flex}) |
           borderHeavy;
  });

  screen.Loop(renderer);
  return 0;
}
