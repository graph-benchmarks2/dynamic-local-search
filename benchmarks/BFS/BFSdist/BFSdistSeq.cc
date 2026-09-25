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

#include "BFSdistSeq.h"

namespace gbbs {

template <class Graph>
double BFSdistSeq_runner(Graph& G, commandLine P) {
  uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
  double dist_cap =
      P.getOptionDoubleValue("-cap", std::numeric_limits<double>::infinity());

  std::cout << "### Application: BFSdistSeq" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: 1" << std::endl;
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
  auto dists = BFSdistSeq(G, src, dist_cap);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

generate_main(gbbs::BFSdistSeq_runner, false);
