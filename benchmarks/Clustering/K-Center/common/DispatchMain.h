// This is a generic, shared dispatch for algorithm benchmark binaries that run on both weighted and unweighted graphs 
// that execute different behavior and subroutines based on the input type.
//
// Typical usage from an algorithm .cc:
//
//   int main(int argc, char* argv[]) {
//     auto app = [](auto& G, gbbs::commandLine P) {
//       return gbbs::MyAlgorithm_runner(G, P);
//     };
//     return gbbs::kcenter_common::dispatch_main(argc, argv, app, false);
//   }

#pragma once

#include <cstddef>
#include <utility>

#include "gbbs/benchmark.h"
#include "gbbs/gbbs.h"
#include "parlay/parallel.h"
#include "benchmarks/Clustering/K-Center/common/GraphFormatGuard.h"

namespace gbbs {
namespace kcenter_common {

template <class App>
int dispatch_main(int argc, char* argv[], App app, bool mutates) {
  gbbs::commandLine P(argc, argv, " [-s] <inFile>");
  bool single_core = P.getOption("-sc");

  auto run = [&]() -> int {
    char* iFile = P.getArgument(0);
    bool symmetric = P.getOptionValue("-s");
    bool compressed = P.getOptionValue("-c");
    bool mmap = P.getOptionValue("-m");
    bool binary = P.getOptionValue("-b");
    size_t rounds = P.getOptionLongValue("-rounds", 3);

    auto format =
        gbbs::kcenter_common::detect_graph_weight_format_or_die(
            iFile, compressed, binary);

    double total_time = 0.0;

    if (format == gbbs::kcenter_common::GraphWeightFormat::kWeighted) {
      if (compressed) {
        if (symmetric) {
          auto G =
              gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::intE>(
                  iFile, mmap);
          for (size_t r = 0; r < rounds; ++r) {
            if (mutates) {
              auto G_copy = G;
              total_time += app(G_copy, P);
            } else {
              total_time += app(G, P);
            }
          }
        } else {
          auto G =
              gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::intE>(
                  iFile, mmap);
          for (size_t r = 0; r < rounds; ++r) {
            if (mutates) {
              auto G_copy = G;
              total_time += app(G_copy, P);
            } else {
              total_time += app(G, P);
            }
          }
        }
      } else {
        if (symmetric) {
          auto G =
              gbbs::gbbs_io::read_weighted_symmetric_graph<gbbs::intE>(
                  iFile, mmap, binary);
          for (size_t r = 0; r < rounds; ++r) {
            if (mutates) {
              auto G_copy = G;
              total_time += app(G_copy, P);
            } else {
              total_time += app(G, P);
            }
          }
        } else {
          auto G =
              gbbs::gbbs_io::read_weighted_asymmetric_graph<gbbs::intE>(
                  iFile, mmap, binary);
          for (size_t r = 0; r < rounds; ++r) {
            if (mutates) {
              auto G_copy = G;
              total_time += app(G_copy, P);
            } else {
              total_time += app(G, P);
            }
          }
        }
      }
    } else {
      if (compressed) {
        if (symmetric) {
          auto G =
              gbbs::gbbs_io::read_compressed_symmetric_graph<gbbs::empty>(
                  iFile, mmap);
          for (size_t r = 0; r < rounds; ++r) {
            if (mutates) {
              auto G_copy = G;
              total_time += app(G_copy, P);
            } else {
              total_time += app(G, P);
            }
          }
        } else {
          auto G =
              gbbs::gbbs_io::read_compressed_asymmetric_graph<gbbs::empty>(
                  iFile, mmap);
          for (size_t r = 0; r < rounds; ++r) {
            if (mutates) {
              auto G_copy = G;
              total_time += app(G_copy, P);
            } else {
              total_time += app(G, P);
            }
          }
        }
      } else {
        if (symmetric) {
          auto G = gbbs::gbbs_io::read_unweighted_symmetric_graph(
              iFile, mmap, binary);
          for (size_t r = 0; r < rounds; ++r) {
            if (mutates) {
              auto G_copy = G;
              total_time += app(G_copy, P);
            } else {
              total_time += app(G, P);
            }
          }
        } else {
          auto G = gbbs::gbbs_io::read_unweighted_asymmetric_graph(
              iFile, mmap, binary);
          for (size_t r = 0; r < rounds; ++r) {
            if (mutates) {
              auto G_copy = G;
              total_time += app(G_copy, P);
            } else {
              total_time += app(G, P);
            }
          }
        }
      }
    }

    auto time_per_iter = total_time / rounds;
    std::cout << "# time per iter: " << time_per_iter << "\n";
    return 0;
  };

  if (!single_core) {
    return run();
  }

  int result = 0;
  parlay::execute_with_scheduler(1, [&]() {
    result = run();
  });
  return result;
}

}  // namespace kcenter_common
}  // namespace gbbs
