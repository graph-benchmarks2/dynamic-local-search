// This file provides an implementation of Gonzalez's (2+epsilon) k-center approximation algorithm.

#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>

#include "gbbs/gbbs.h"
#include "benchmarks/BFS/BFSdist/BFSdist.h"
#include "benchmarks/BFS/BFSdist/BFSdistRepair.h"
#include "benchmarks/BFS/BFSdist/BFSdistSeq.h"
#include "benchmarks/BFS/BFSdist/BFSdistRepairSeq.h"
#include "benchmarks/Clustering/K-Center/Gonzalez/common/KCenterResult.h"
#include "benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaStepping.h"
#include "benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaSteppingRepair.h"
#include "benchmarks/PositiveWeightSSSP/Dijkstra/Dijkstra.h"
#include "benchmarks/PositiveWeightSSSP/Dijkstra/DijkstraRepair.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "parlay/sequence.h"

namespace gbbs {

template <class Graph>
KCenterResult<float> Gonzalez(
    Graph& G,
    uintE k,
    double delta,
    size_t num_buckets = 128,
    uint64_t seed = 12345,
    bool single_core = false,
    bool sssp_verbose = false) {

  using W = typename Graph::weight_type;
  using vertex_id = uintE;
  using Distance = float;

  const size_t n = G.n;
  const Distance kInf = std::numeric_limits<Distance>::max();

  parlay::sequence<Distance> dist(n, kInf);
  parlay::sequence<vertex_id> centers;

  if (k == 0 || n == 0) {
    return {std::move(centers), static_cast<Distance>(0), n};
  }

  // Pick initial center uniformly at random
  parlay::random rng(seed);
  uint64_t r = rng.ith_rand(0);
  vertex_id first = static_cast<vertex_id>(r % n);
  centers.push_back(first);

  // First center: run appropriate from-scratch SSSP
  if constexpr (std::is_same<W, gbbs::empty>::value) {
    auto initial_distances =
    single_core
        ? BFSdistSeq(
              G,
              first,
              std::numeric_limits<double>::infinity(),
              sssp_verbose)
        : BFSdist(
              G,
              first,
              std::numeric_limits<double>::infinity(),
              sssp_verbose);
    parallel_for(0, n, [&](size_t i) {
      dist[i] = (initial_distances[i] == UINT_E_MAX)
                    ? kInf
                    : static_cast<Distance>(initial_distances[i]);
    });
  } else {
    auto initial_distances =
    single_core
        ? Dijkstra(
              G,
              first,
              std::numeric_limits<double>::infinity(),
              sssp_verbose)
        : DeltaStepping(
              G,
              first,
              delta,
              num_buckets,
              std::numeric_limits<double>::infinity(),
              sssp_verbose);

    parallel_for(0, n, [&](size_t i) {
      dist[i] = static_cast<Distance>(initial_distances[i]);
    });
  }

  // Remaining centers use previously computed distances
  for (size_t i = 1; i < k; i++) {
    auto max_pair_monoid = parlay::make_monoid(
        [](auto a, auto b) { return (a.second > b.second) ? a : b; },
        std::make_pair(size_t{0}, std::numeric_limits<Distance>::lowest()));

    auto max_pair = parlay::reduce(
        parlay::tabulate(n, [&](size_t j) {
          return std::make_pair(j, dist[j]);
        }),
        max_pair_monoid);

    auto [farthest_vtx, max_dist] = max_pair;
    (void)max_dist;

    centers.push_back(static_cast<vertex_id>(farthest_vtx));

    if constexpr (std::is_same<W, gbbs::empty>::value) {
      auto seeds = parlay::sequence<std::pair<uintE, uintE>>(1);
      seeds[0] = std::make_pair(static_cast<uintE>(farthest_vtx),
                                static_cast<uintE>(0));

      auto dist_uint = parlay::sequence<uintE>(n, UINT_E_MAX);
      parallel_for(0, n, [&](size_t j) {
        dist_uint[j] = (dist[j] == kInf) ? UINT_E_MAX
                                         : static_cast<uintE>(dist[j]);
      });

      if (single_core) {
        BFSdistRepairSeq(G, dist_uint, seeds);
      } else {
        BFSdistRepair(G, dist_uint, seeds);
      }

      parallel_for(0, n, [&](size_t j) {
        dist[j] = (dist_uint[j] == UINT_E_MAX)
                      ? kInf
                      : static_cast<Distance>(dist_uint[j]);
      });
    } else {
      auto seeds = parlay::sequence<std::pair<uintE, Distance>>(1);
      seeds[0] =
          std::make_pair(static_cast<uintE>(farthest_vtx),
                         static_cast<Distance>(0));

      if (single_core) {
        DijkstraRepair(G, dist, seeds);
      } else {
        DeltaSteppingRepair(G, dist, seeds, delta, num_buckets);
      }
    }
  }

  size_t unreachable_vertices = 0;
  Distance max_dist_to_centers = static_cast<Distance>(0);

  for (size_t i = 0; i < n; ++i) {
    if (dist[i] == kInf) {
      unreachable_vertices++;
    } else {
      max_dist_to_centers = std::max(max_dist_to_centers, dist[i]);
    }
  }

  return {std::move(centers), max_dist_to_centers, unreachable_vertices};
}

}  // namespace gbbs