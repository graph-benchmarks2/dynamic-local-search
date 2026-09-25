// Exact insertion only APSP maintenance for undirected graphs. The batch of distance decreases is fully processed and the metric is restored before being passed on to the caller.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "gbbs/gbbs.h"

#include "benchmarks/BFS/BFSdist/BFSdistRepairSeq.h"
#include "benchmarks/BFS/BFSdist/BFSdistSeq.h"
#include "benchmarks/Clustering/K-Median/common/DistanceDecrease.h"
#include "benchmarks/Clustering/K-Median/common/KMedianCommon.h"
#include "benchmarks/PositiveWeightSSSP/Dijkstra/Dijkstra.h"
#include "benchmarks/PositiveWeightSSSP/Dijkstra/DijkstraRepair.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"

namespace gbbs {
namespace kmedian {

template <class Graph>
class IncrementalAPSP {
 public:
  using W = typename Graph::weight_type;
  using Distance =
      typename std::conditional<std::is_same<W, gbbs::empty>::value,
                                uintE, W>::type;

  static constexpr Distance kMaxDistance =
      std::numeric_limits<Distance>::max();

  IncrementalAPSP() = default;

  explicit IncrementalAPSP(
      Graph& G,
      double dist_cap = std::numeric_limits<double>::infinity(),
      SSSPTiming* timing = nullptr) {
    Build(G, dist_cap, timing);
  }

  void Build(
      Graph& G,
      double dist_cap = std::numeric_limits<double>::infinity(),
      SSSPTiming* timing = nullptr) {
    n_ = G.n;
    dist_cap_ = dist_cap;

    if (n_ == 0) {
      Fail("IncrementalAPSP requires a nonempty graph.");
    }

    distances_.clear();
    distances_.resize(n_);

    std::vector<double> source_times(n_, 0.0);

    parlay::parallel_for(0, n_, [&](size_t source) {
      timer sssp_timer;
      sssp_timer.start();

      if constexpr (std::is_same<W, gbbs::empty>::value) {
        distances_[source] =
            BFSdistSeq(
                G,
                static_cast<uintE>(source),
                dist_cap_,
                /*verbose=*/false);
      } else {
        distances_[source] =
            Dijkstra(
                G,
                static_cast<uintE>(source),
                dist_cap_,
                /*verbose=*/false);
      }

      source_times[source] = sssp_timer.stop();
    });

    if (timing != nullptr) {
      for (double elapsed : source_times) {
        timing->Add(elapsed);
      }
    }

    bfs_scratch_ = BFSdistRepairSeqScratch{};
    dijkstra_scratch_ = DijkstraRepairScratch{};
  }

  size_t n() const {
    return n_;
  }

  Distance distance(uintE u, uintE v) const {
    CheckVertex(u);
    CheckVertex(v);
    return distances_[u][v];
  }

  void set_distance(uintE u, uintE v, double value) {
    CheckVertex(u);
    CheckVertex(v);

    const Distance cast_value =
        static_cast<Distance>(value);

    distances_[u][v] = cast_value;
    distances_[v][u] = cast_value;
  }

  const sequence<Distance>& row(uintE source) const {
    CheckVertex(source);
    return distances_[source];
  }

  parlay::sequence<DistanceDecrease> GenerateBatch(
      Graph& G,
      uintE u,
      uintE v,
      W w = W{}) {
    CheckGraph(G);
    CheckVertex(u);
    CheckVertex(v);

    if (u == v) {
      Fail("IncrementalAPSP does not support self-loop insertions.");
    }

    const Distance ww = EdgeWeight(w);

    if constexpr (!std::is_same<W, gbbs::empty>::value) {
      if (!(ww > static_cast<Distance>(0))) {
        Fail("IncrementalAPSP requires positive inserted edge weights.");
      }
    }

    struct RollbackEntry {
      uintE a;
      uintE b;
      Distance old_distance;
      Distance new_distance;
    };

    std::vector<RollbackEntry> changed_pairs;
    std::vector<uintE> affected;

    for (uintE source = 0; source < n_; ++source) {
      sequence<std::pair<uintE, Distance>> seeds;

      const Distance du = distances_[source][u];
      const Distance dv = distances_[source][v];

      const Distance via_uv = SafeAdd(du, ww);
      const Distance via_vu = SafeAdd(dv, ww);

      const bool improve_v = via_uv < dv;
      const bool improve_u = via_vu < du;

      if (!improve_v && !improve_u) {
        continue;
      }

      seeds.reserve(
          static_cast<size_t>(improve_v) +
          static_cast<size_t>(improve_u));

      if (improve_v) {
        seeds.push_back(std::make_pair(v, via_uv));
      }

      if (improve_u) {
        seeds.push_back(std::make_pair(u, via_vu));
      }

      if constexpr (std::is_same<W, gbbs::empty>::value) {
        BFSdistRepairSeq(
            G,
            distances_[source],
            seeds,
            dist_cap_,
            nullptr,
            affected,
            bfs_scratch_);
      } else {
        DijkstraRepair(
            G,
            distances_[source],
            seeds,
            dist_cap_,
            static_cast<AlphaBucketState<Distance>*>(nullptr),
            affected,
            dijkstra_scratch_);
      }

      for (uintE target : affected) {
        if (source >= target) {
          continue;
        }

        const Distance new_distance =
            distances_[source][target];

        const Distance old_distance =
            distances_[target][source];

        if (!(new_distance < old_distance)) {
          Fail("IncrementalAPSP detected an asymmetric or non-decreasing update.");
        }

        changed_pairs.push_back(
            RollbackEntry{
                source,
                target,
                old_distance,
                new_distance});
      }
    }

    for (const auto& change : changed_pairs) {
      distances_[change.a][change.b] =
          change.old_distance;

      distances_[change.b][change.a] =
          change.old_distance;
    }

    std::sort(
        changed_pairs.begin(),
        changed_pairs.end(),
        [](const RollbackEntry& lhs, const RollbackEntry& rhs) {
          if (lhs.a != rhs.a) return lhs.a < rhs.a;
          return lhs.b < rhs.b;
        });

    parlay::sequence<DistanceDecrease> batch(
        changed_pairs.size());

    for (size_t i = 0; i < changed_pairs.size(); ++i) {
      batch[i] =
          DistanceDecrease{
              changed_pairs[i].a,
              changed_pairs[i].b,
              static_cast<double>(
                  changed_pairs[i].new_distance)};
    }

    return batch;
  }

 private:
  Distance EdgeWeight(const W& w) const {
    if constexpr (std::is_same<W, gbbs::empty>::value) {
      return static_cast<Distance>(1);
    } else {
      return static_cast<Distance>(w);
    }
  }

  Distance SafeAdd(
      Distance d,
      Distance w) const {
    if (d == kMaxDistance) {
      return kMaxDistance;
    }

    if (d > kMaxDistance - w) {
      return kMaxDistance;
    }

    return d + w;
  }

  void CheckVertex(uintE v) const {
    if (v >= n_) {
      Fail("IncrementalAPSP vertex is outside the valid range.");
    }
  }

  void CheckGraph(const Graph& G) const {
    if (G.n != n_) {
      Fail("IncrementalAPSP graph size changed.");
    }
  }

  [[noreturn]] static void Fail(const char* message) {
    std::cout << "ERROR: " << message << std::endl;
    std::exit(-1);
  }

  size_t n_ = 0;
  double dist_cap_ =
      std::numeric_limits<double>::infinity();

  std::vector<sequence<Distance>> distances_;

  BFSdistRepairSeqScratch bfs_scratch_;
  DijkstraRepairScratch dijkstra_scratch_;
};

}  // namespace kmedian
}  // namespace gbbs
