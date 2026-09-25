// flags:
//   required:
//     -src <vertex> : source vertex
//   optional:
//     -rounds <int> : number of times to run the algorithm
//     -s            : input graph is symmetric / undirected
//     -c            : input graph is compressed
//     -m            : input graph should be mmap'd
//     -b            : input graph is binary
//     -cap <dbl>    : distance cap; relaxations above this value are ignored

#define WEIGHTED 1

#include "DijkstraRepair.h"

namespace gbbs {

template <class Graph>
double DijkstraRepair_runner(Graph& G, commandLine P) {
  using W = typename Graph::weight_type;
  using Distance =
      typename std::conditional<std::is_same<W, gbbs::empty>::value, uintE,
                                W>::type;

  uintE src = P.getOptionLongValue("-src", 0);
  double dist_cap =
      P.getOptionDoubleValue("-cap", std::numeric_limits<double>::infinity());

  std::cout << "### Application: DijkstraRepair" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: 1" << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -src = " << src
            << " -cap = " << dist_cap << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  constexpr Distance kMaxWeight = std::numeric_limits<Distance>::max();
  auto dist = sequence<Distance>(G.n, kMaxWeight);

  auto seeds = sequence<std::pair<uintE, Distance>>(1);
  seeds[0] = std::make_pair(src, static_cast<Distance>(0));

  timer t;
  t.start();
  DijkstraRepair(G, dist, seeds, dist_cap);
  double tt = t.stop();

  auto dist_im = parlay::delayed_seq<Distance>(
      G.n, [&](size_t i) { return (dist[i] == kMaxWeight) ? 0 : dist[i]; });

  std::cout << "max_dist = " << parlay::reduce_max(dist_im) << std::endl;
  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

generate_weighted_main(gbbs::DijkstraRepair_runner, false);