// flags:
//   required:
//     -src <vertex> : initial source vertex
//   optional:
//     -rounds <int>       : number of times to run the algorithm
//     -s                  : input graph is symmetric / undirected
//     -m                  : input graph should be mmap'd
//     -c                  : compressed graph input
//     -b                  : binary graph input
//     -sc                 : use single-core SSSP backend
//     -delta <dbl>        : bucket width 
//     -nb <int>           : number of buckets; must be a power of two
//     -cap <dbl>          : distance cap
//     -center <vertex>    : after initialization, insert this vertex as an
//                           additional zero-distance center
//     -u <vertex>         : source endpoint of an inserted edge
//     -v <vertex>         : target endpoint of an inserted edge
//     -w <weight>         : weight of inserted edge for weighted graphs;
//                           ignored for unweighted graphs, where weight is 1
//     -insert_undirected  : insert both (u,v) and (v,u)

#include <algorithm>
#include <limits>
#include <type_traits>

#include "gbbs/gbbs.h"
#include "benchmarks/Clustering/K-Center/common/DispatchMain.h"
#include "RR.h"

namespace gbbs {

template <class Graph>
double RR_runner(Graph& G, commandLine P) {
  using W = typename Graph::weight_type;
  using Distance = typename RR<Graph>::Distance;

  uintE src = P.getOptionLongValue("-src", 0);
  size_t num_buckets = P.getOptionLongValue("-nb", 32);
  double delta = P.getOptionDoubleValue("-delta", 1.0);
  double dist_cap =
      P.getOptionDoubleValue("-cap", std::numeric_limits<double>::infinity());
  bool single_core = P.getOption("-sc");

  long center_opt = P.getOptionLongValue("-center", -1);
  long u_opt = P.getOptionLongValue("-u", -1);
  long v_opt = P.getOptionLongValue("-v", -1);
  long w_opt = P.getOptionLongValue("-w", -1);
  bool undirected = P.getOption("-insert_undirected");

  std::cout << "### Application: RR" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: "
          << (single_core ? 1 : num_workers())
          << std::endl;
  std::cout << "### Graph Type: "
            << (std::is_same<W, gbbs::empty>::value ? "unweighted" : "weighted")
            << std::endl;
  std::cout << "### SSSP Mode: "
            << (single_core ? "single-core" : "parallel") << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -src = " << src
            << ", -delta = " << delta
            << ", -nb = " << num_buckets
            << ", -cap = " << dist_cap
            << (single_core ? ", -sc = true" : ", -sc = false")
            << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  if constexpr (!std::is_same<W, gbbs::empty>::value) {
    if (!single_core &&
        num_buckets != (((uintE)1) << parlay::log2_up(num_buckets))) {
      std::cout << "Please specify a number of buckets that is a power of two\n";
      exit(-1);
    }
  }

  RR<Graph> rr(G, src, delta, num_buckets, dist_cap, single_core);

  timer t_init;
  t_init.start();
  rr.initialize();
  double init_time = t_init.stop();

  auto dist0 = rr.distances();
  auto dist0_im = parlay::delayed_seq<Distance>(
      G.n, [&](size_t i) {
        return (dist0[i] == RR<Graph>::kMaxDistance) ? 0 : dist0[i];
      });

  std::cout << "initial max_dist = "
            << parlay::reduce_max(dist0_im) << std::endl;
  std::cout << "initial support edges = "
            << rr.support_out_edges().size() << std::endl;
  std::cout << "### Init Time: " << init_time << std::endl;

  if (center_opt >= 0) {
    uintE c = static_cast<uintE>(center_opt);

    std::cout << "### Update: insert_center(" << c << ")" << std::endl;

    timer t_upd;
    t_upd.start();
    bool did_update = rr.insert_center(c);
    double upd_time = t_upd.stop();

    auto dist1 = rr.distances();
    auto dist1_im = parlay::delayed_seq<Distance>(
        G.n, [&](size_t i) {
          return (dist1[i] == RR<Graph>::kMaxDistance) ? 0 : dist1[i];
        });

    std::cout << "update changed distances = " << did_update << std::endl;
    std::cout << "updated max_dist = "
              << parlay::reduce_max(dist1_im) << std::endl;
    std::cout << "updated support edges = "
              << rr.support_out_edges().size() << std::endl;
    std::cout << "### Update Time: " << upd_time << std::endl;
    return upd_time;
  }

  if (u_opt >= 0 || v_opt >= 0 || w_opt >= 0) {
    if (u_opt < 0 || v_opt < 0) {
      std::cout << "To insert an edge, please specify -u and -v." << std::endl;
      exit(-1);
    }

    uintE u = static_cast<uintE>(u_opt);
    uintE v = static_cast<uintE>(v_opt);

    W w{};
    if constexpr (std::is_same<W, gbbs::empty>::value) {
      if (w_opt >= 0) {
        std::cout << "### Note: ignoring -w for unweighted graph; edge weight is 1."
                  << std::endl;
      }
    } else {
      if (w_opt < 0) {
        std::cout << "To insert a weighted edge, please specify -w." << std::endl;
        exit(-1);
      }
      w = static_cast<W>(w_opt);
    }

    std::cout << "### Update: "
              << (undirected ? "insert_undirected_edge(" : "insert_edge(")
              << u << ", " << v;
    if constexpr (!std::is_same<W, gbbs::empty>::value) {
      std::cout << ", " << w_opt;
    }
    std::cout << ")" << std::endl;

    timer t_upd;
    t_upd.start();
    bool did_update = undirected ? rr.insert_undirected_edge(u, v, w)
                                 : rr.insert_edge(u, v, w);
    double upd_time = t_upd.stop();

    auto dist1 = rr.distances();
    auto dist1_im = parlay::delayed_seq<Distance>(
        G.n, [&](size_t i) {
          return (dist1[i] == RR<Graph>::kMaxDistance) ? 0 : dist1[i];
        });

    std::cout << "update changed distances = " << did_update << std::endl;
    std::cout << "updated max_dist = "
              << parlay::reduce_max(dist1_im) << std::endl;
    std::cout << "updated support edges = "
              << rr.support_out_edges().size() << std::endl;
    std::cout << "### Update Time: " << upd_time << std::endl;
    return upd_time;
  }

  return init_time;
}

}  // namespace gbbs

int main(int argc, char* argv[]) {
  auto app = [](auto& G, gbbs::commandLine P) {
    return gbbs::RR_runner(G, P);
  };
  return gbbs::kcenter_common::dispatch_main(argc, argv, app, false);
}