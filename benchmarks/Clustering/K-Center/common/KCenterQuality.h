#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "gbbs/gbbs.h"
#include "benchmarks/BFS/BFSdist/BFSdist.h"
#include "benchmarks/BFS/BFSdist/BFSdistSeq.h"
#include "benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaStepping.h"
#include "benchmarks/PositiveWeightSSSP/Dijkstra/Dijkstra.h"
#include "parlay/sequence.h"

namespace gbbs {
namespace kcenter_quality {

template <class Graph>
using KCenterDistance =
    typename std::conditional<std::is_same<typename Graph::weight_type,
                                           gbbs::empty>::value,
                              uintE,
                              typename Graph::weight_type>::type;

template <class Distance>
struct Result {
  Distance max_finite_dist_to_centers = static_cast<Distance>(0);
  size_t unreachable_vertices = 0;

  bool has_unreachable() const {
    return unreachable_vertices > 0;
  }

  std::string max_dist_to_centers_string() const {
    if (has_unreachable()) return "inf";
    return std::to_string(max_finite_dist_to_centers);
  }
};

template <class Distance>
inline Distance infinity() {
  return std::numeric_limits<Distance>::max();
}

template <class RawDistance>
inline bool is_unreachable_raw(RawDistance x) {
  using R = typename std::decay<RawDistance>::type;

  if constexpr (std::is_floating_point<R>::value) {
    if (std::isinf(x)) return true;

    // Some GBBS weighted SSSP paths use/cast uintE max as the unreachable
    // sentinel, which appears as about 2.14748e+09 for float.
    if (x == static_cast<R>(std::numeric_limits<uintE>::max())) return true;

    return x == std::numeric_limits<R>::max();
  } else {
    return x == std::numeric_limits<R>::max();
  }
}

template <class Graph>
parlay::sequence<typename Graph::edge> extract_unique_undirected_edges(Graph& G) {
  using W = typename Graph::weight_type;
  using edge = typename Graph::edge;

  std::vector<edge> edges;
  edges.reserve(G.m / 2 + 1);

  for (uintE u = 0; u < G.n; ++u) {
    auto map_f = [&](const uintE& src, const uintE& v, const W& w) {
      if (src <= v) edges.emplace_back(src, v, w);
    };
    G.get_vertex(u).out_neighbors().map(map_f, false);
  }

  return parlay::sequence<edge>(edges.begin(), edges.end());
}

template <class W>
using TempSymGraph = gbbs::symmetric_graph<gbbs::symmetric_vertex, W>;

template <class Graph>
TempSymGraph<typename Graph::weight_type> build_augmented_graph_with_supersource(
    const parlay::sequence<typename Graph::edge>& base_edges,
    size_t n,
    const parlay::sequence<uintE>& centers) {
  using W = typename Graph::weight_type;
  using edge = std::tuple<uintE, uintE, W>;

  const uintE s_star = static_cast<uintE>(n);
  parlay::sequence<edge> edges_aug(base_edges.size() + centers.size());

  parallel_for(0, base_edges.size(), [&](size_t i) {
    edges_aug[i] = base_edges[i];
  });

  parallel_for(0, centers.size(), [&](size_t i) {
    edges_aug[base_edges.size() + i] =
        std::make_tuple(s_star, centers[i], W{});
  });

  return TempSymGraph<W>::from_edges(edges_aug, n + 1);
}

template <class Graph>
Result<KCenterDistance<Graph>> compute(
    Graph& G,
    const parlay::sequence<uintE>& centers,
    double delta,
    size_t num_buckets,
    bool single_core = false) {
  using W = typename Graph::weight_type;
  using Distance = KCenterDistance<Graph>;

  Result<Distance> out;

  if (G.n == 0) return out;

  if (centers.empty()) {
    out.unreachable_vertices = G.n;
    return out;
  }

  auto base_edges = extract_unique_undirected_edges(G);
  auto G_aug = build_augmented_graph_with_supersource<Graph>(
      base_edges, G.n, centers);
  uintE s_star = static_cast<uintE>(G.n);

  if constexpr (std::is_same<W, gbbs::empty>::value) {
    auto dist_aug =
        single_core
            ? BFSdistSeq(G_aug, s_star,
                         std::numeric_limits<double>::infinity(),
                         /*verbose=*/false)
            : BFSdist(G_aug, s_star,
                      std::numeric_limits<double>::infinity(),
                      /*verbose=*/false);

    for (size_t i = 0; i < G.n; ++i) {
      auto raw = dist_aug[i];

      if (is_unreachable_raw(raw)) {
        out.unreachable_vertices++;
      } else {
        Distance d = raw == 0
                         ? static_cast<Distance>(0)
                         : static_cast<Distance>(raw - 1);
        out.max_finite_dist_to_centers =
            std::max(out.max_finite_dist_to_centers, d);
      }
    }
  } else {
    auto dist_aug =
        single_core
            ? Dijkstra(G_aug, s_star,
                       std::numeric_limits<double>::infinity(),
                       /*verbose=*/false)
            : DeltaStepping(G_aug, s_star,
                            delta,
                            num_buckets,
                            std::numeric_limits<double>::infinity(),
                            /*verbose=*/false);


    size_t raw_inf = 0;
    size_t raw_uintE_inf = 0;
    size_t raw_big = 0;

    for (size_t i = 0; i < G.n; ++i) {
    auto raw = dist_aug[i];

    if (raw == std::numeric_limits<decltype(raw)>::max()) raw_inf++;
    if (raw == static_cast<decltype(raw)>(std::numeric_limits<uintE>::max())) raw_uintE_inf++;
    if (static_cast<double>(raw) > 1e8) raw_big++;
    }

    //std::cout << "DEBUG quality raw_inf = " << raw_inf << std::endl;
    //std::cout << "DEBUG quality raw_uintE_inf = " << raw_uintE_inf << std::endl;
    //std::cout << "DEBUG quality raw_big = " << raw_big << std::endl;  

    for (size_t i = 0; i < G.n; ++i) {
      auto raw = dist_aug[i];

      if (is_unreachable_raw(raw)) {
        out.unreachable_vertices++;
      } else {
        Distance d = static_cast<Distance>(raw);
        out.max_finite_dist_to_centers =
            std::max(out.max_finite_dist_to_centers, d);
      }
    }
  }

  return out;
}

template <class Distance>
void print(const Result<Distance>& q) {
  std::cout << "max_finite_dist_to_centers = "
            << q.max_finite_dist_to_centers << std::endl;
  std::cout << "max_dist_to_centers = "
            << q.max_dist_to_centers_string() << std::endl;
  std::cout << "unreachable_vertices = "
            << q.unreachable_vertices << std::endl;
}

}  // namespace kcenter_quality
}  // namespace gbbs