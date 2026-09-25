// flags:
//   required:
//     -src <vertex> : source vertex used as initial repair seed
//   optional:
//     -rounds <int> : number of times to run the algorithm
//     -s            : input graph is symmetric / undirected
//     -c            : input graph is compressed
//     -m            : input graph should be mmap'd
//     -b            : input graph is binary
//     -cap <dbl>    : hop-distance cap; relaxations above this value are ignored

#include "BFSdistRepair.h"

namespace gbbs {

template <class Graph>
double BFSdistRepair_runner(Graph& G, commandLine P) {
  uintE src = static_cast<uintE>(P.getOptionLongValue("-src", 0));
  double dist_cap =
      P.getOptionDoubleValue("-cap", std::numeric_limits<double>::infinity());

  std::cout << "### Application: BFSdistRepair" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -src = " << src
            << " -cap = " << dist_cap << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  auto dist = sequence<uintE>(G.n, UINT_E_MAX);

  auto seeds = sequence<std::pair<uintE, uintE>>(1);
  seeds[0] = std::make_pair(src, static_cast<uintE>(0));

  timer t;
  t.start();
  BFSdistRepair(G, dist, seeds, dist_cap);
  double tt = t.stop();

  auto dist_im = parlay::delayed_seq<uintE>(
      G.n, [&](size_t i) { return (dist[i] == UINT_E_MAX) ? 0 : dist[i]; });

  std::cout << "max_dist = " << parlay::reduce_max(dist_im) << std::endl;
  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

generate_main(gbbs::BFSdistRepair_runner, false);