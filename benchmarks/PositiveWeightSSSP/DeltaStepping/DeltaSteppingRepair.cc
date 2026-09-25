// flags:
//   required:
//     -src <vertex> : source vertex
//   optional:
//     -rounds <int> : number of times to run the algorithm
//     -s            : input graph is symmetric / undirected
//     -c            : input graph is compressed
//     -m            : input graph should be mmap'd
//     -b            : input graph is binary
//     -delta <dbl>  : bucket width
//     -nb <int>     : number of buckets; must be a power of two
//     -cap <dbl>    : distance cap; relaxations above this value are ignored

#define WEIGHTED 1

#include "DeltaSteppingRepair.h"

namespace gbbs {

template <class Graph>
double DeltaSteppingRepair_runner(Graph& G, commandLine P) {
  using W = typename Graph::weight_type;
  using Distance =
      typename std::conditional<std::is_same<W, gbbs::empty>::value, uintE, W>::type;

  uintE src = P.getOptionLongValue("-src", 0);
  size_t num_buckets = P.getOptionLongValue("-nb", 32);
  double delta = P.getOptionDoubleValue("-delta", 1.0);
  double dist_cap =
      P.getOptionDoubleValue("-cap", std::numeric_limits<double>::infinity());

  std::cout << "### Application: DeltaSteppingRepair" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -src = " << src << " -delta = " << delta
            << " -nb (num_buckets) = " << num_buckets
            << " -cap = " << dist_cap << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  if (num_buckets != (((uintE)1) << parlay::log2_up(num_buckets))) {
    std::cout << "Please specify a number of buckets that is a power of two"
              << "\n";
    exit(-1);
  }

  constexpr Distance kMaxWeight = std::numeric_limits<Distance>::max();
  auto dist = sequence<Distance>(G.n, kMaxWeight);

  auto seeds = sequence<std::pair<uintE, Distance>>(1);
  seeds[0] = std::make_pair(src, static_cast<Distance>(0));

  timer t;
  t.start();
  DeltaSteppingRepair(G, dist, seeds, delta, num_buckets, dist_cap);
  double tt = t.stop();

  auto dist_im = parlay::delayed_seq<Distance>(
      G.n, [&](size_t i) { return (dist[i] == kMaxWeight) ? 0 : dist[i]; });

  std::cout << "max_dist = " << parlay::reduce_max(dist_im) << std::endl;
  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}

}  // namespace gbbs

generate_weighted_main(gbbs::DeltaSteppingRepair_runner, false);