// flags:
//   required:
//     -src <vertex> : source vertex
//   optional:
//     -rounds <int> : number of times to run the algorithm
//     -s            : input graph is symmetric / undirected
//     -c            : input graph is compressed
//     -m            : input graph should be mmap'd
//     -b            : input graph is binary
//     -cap <dbl>    : hop-distance cap; relaxations above this value are ignored

#include "BFSdist.h"

namespace gbbs {

template <class Graph>
double BFSdist_runner(Graph& G, commandLine P) {
  uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
  double dist_cap =
      P.getOptionDoubleValue("-cap", std::numeric_limits<double>::infinity());

  std::cout << "### Application: BFSdist" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -src = " << src;
  if (!std::isinf(dist_cap)) {
    std::cout << " -cap = " << dist_cap;
  }
  std::cout << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  timer t;
  t.start();
  auto dists = BFSdist(G, src, dist_cap);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

generate_main(gbbs::BFSdist_runner, false);