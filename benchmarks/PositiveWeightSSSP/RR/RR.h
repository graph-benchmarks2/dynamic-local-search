// This file provides an implementation of the dynamic SSSP algorithm by Ramalingam and Reps (incremental only at the moment)

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

#include "gbbs/gbbs.h"
#include "gbbs/vertex_subset.h"
#include "gbbs/edge_map_data.h"

#include "benchmarks/BFS/BFSdist/BFSdist.h"
#include "benchmarks/BFS/BFSdist/BFSdistRepair.h"
#include "benchmarks/BFS/BFSdist/BFSdistSeq.h"
#include "benchmarks/BFS/BFSdist/BFSdistRepairSeq.h"

#include "benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaStepping.h"
#include "benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaSteppingRepair.h"
#include "benchmarks/PositiveWeightSSSP/Dijkstra/Dijkstra.h"
#include "benchmarks/PositiveWeightSSSP/Dijkstra/DijkstraRepair.h"

namespace gbbs {

// Different strategies to build support graph structure:  
// - Rebuild->wait for SSSP update to finish and build entire SPG after
// - Local->Mark affected vertices during SSSP and modify SPG after
// - Hybrid->Start off with Rebuild mode and switch to local at a certain threshold
enum class RRSPGMode {
  Rebuild,
  Local,
  Hybrid
};

template <class Graph>
struct RR {
  using W = typename Graph::weight_type;
  using Distance =
      typename std::conditional<std::is_same<W, gbbs::empty>::value,
                                uintE, W>::type;

  static constexpr Distance kMaxDistance =
      std::numeric_limits<Distance>::max();

  Graph& G;
  uintE source;
  double delta;
  size_t num_buckets;
  Distance dist_cap;
  bool single_core;

  RRSPGMode spg_mode;
  double spg_threshold;

  size_t spg_rebuild_calls = 0;
  size_t spg_local_calls = 0;
  size_t spg_hybrid_rebuild_calls = 0;
  size_t spg_hybrid_local_calls = 0;

  size_t spg_total_affected = 0;
  size_t spg_max_affected = 0;

  sequence<Distance> dist;

  struct OverlayEdge {
    uintE to;
    Distance w;
  };

  std::vector<std::vector<OverlayEdge>> overlay_out;

  std::vector<bool> virtual_source_attached;
  sequence<uintE> attached_centers;

  // CSR representation of the support graph.
  sequence<uintE> support_offsets;
  sequence<uintE> support_edges;

  // Mutable adjacency representation used to update the support graph locally.
  std::vector<std::vector<uintE>> support_adj;

  RR(Graph& _G,
     uintE _src,
     double _delta,
     size_t _num_buckets,
     double _dist_cap = std::numeric_limits<double>::infinity(),
     bool _single_core = false,
     RRSPGMode _spg_mode = RRSPGMode::Rebuild,
     double _spg_threshold = 0.25)
      : G(_G),
        source(_src),
        delta(_delta),
        num_buckets(_num_buckets),
        dist_cap((std::isinf(_dist_cap)
                      ? kMaxDistance
                      : static_cast<Distance>(_dist_cap))),
        single_core(_single_core),
        spg_mode(_spg_mode),
        spg_threshold(_spg_threshold) {}

  inline Distance edge_weight(const W& w) const {
    if constexpr (std::is_same<W, gbbs::empty>::value) {
      return static_cast<Distance>(1);
    } else {
      return static_cast<Distance>(w);
    }
  }

  void initialize() {
    if constexpr (std::is_same<W, gbbs::empty>::value) {
      dist = single_core ? BFSdistSeq(G, source, dist_cap)
                         : BFSdist(G, source, dist_cap);
    } else {
      dist = single_core
                 ? Dijkstra(G, source, dist_cap)
                 : DeltaStepping(G, source, delta, num_buckets, dist_cap);
    }

    overlay_out = std::vector<std::vector<OverlayEdge>>(G.n);
    virtual_source_attached = std::vector<bool>(G.n, false);
    attached_centers = sequence<uintE>();
    support_adj = std::vector<std::vector<uintE>>(G.n);

    rebuild_support_graph();
  }

  void initialize_empty() {
    dist = sequence<Distance>(G.n, kMaxDistance);

    overlay_out = std::vector<std::vector<OverlayEdge>>(G.n);
    virtual_source_attached = std::vector<bool>(G.n, false);
    attached_centers = sequence<uintE>();
    support_adj = std::vector<std::vector<uintE>>(G.n);

    rebuild_support_graph();
  }

  const sequence<Distance>& distances() const { return dist; }

  const sequence<uintE>& support_out_edges() const {
    return support_edges;
  }

  const sequence<uintE>& centers() const {
    return attached_centers;
  }

  size_t num_spg_rebuild_calls() const {
    return spg_rebuild_calls;
  }

  size_t num_spg_local_calls() const {
    return spg_local_calls;
  }

  size_t num_spg_hybrid_rebuild_calls() const {
    return spg_hybrid_rebuild_calls;
  }

  size_t num_spg_hybrid_local_calls() const {
    return spg_hybrid_local_calls;
  }

  size_t total_affected_vertices() const {
    return spg_total_affected;
  }

  size_t max_affected_vertices() const {
    return spg_max_affected;
  }

  double avg_affected_vertices() const {
    const size_t calls = spg_rebuild_calls + spg_local_calls;
    if (calls == 0) return 0.0;
    return static_cast<double>(spg_total_affected) /
          static_cast<double>(calls);
  }

  sequence<uintE> repair_from_seeds(
      const sequence<std::pair<uintE, Distance>>& seeds) {
    if constexpr (std::is_same<W, gbbs::empty>::value) {
      if (single_core) {
        return BFSdistRepairSeqAffected(G, dist, seeds, dist_cap);
      } else {
        return BFSdistRepairAffected(G, dist, seeds, dist_cap);
      }
    } else {
      if (single_core) {
        return DijkstraRepairAffected(G, dist, seeds, dist_cap);
      } else {
        return DeltaSteppingRepairAffected(G, dist, seeds, delta, num_buckets,
                                           dist_cap);
      }
    }
  }

  void maintain_support_graph_after_repair(
      const sequence<uintE>& affected_vertices) {
    if (affected_vertices.empty()) return;

    const size_t affected = affected_vertices.size();
    spg_total_affected += affected;
    spg_max_affected = std::max(spg_max_affected, affected);

    if (spg_mode == RRSPGMode::Rebuild) {
      ++spg_rebuild_calls;
      rebuild_support_graph();
      return;
    }

    if (spg_mode == RRSPGMode::Local) {
      ++spg_local_calls;
      update_support_graph_around(affected_vertices);
      return;
    }

    const double affected_fraction =
        static_cast<double>(affected) /
        static_cast<double>(std::max<size_t>(1, G.n));

    if (affected_fraction > spg_threshold) {
      ++spg_rebuild_calls;
      ++spg_hybrid_rebuild_calls;
      rebuild_support_graph();
    } else {
      ++spg_local_calls;
      ++spg_hybrid_local_calls;
      update_support_graph_around(affected_vertices);
    }
  }

  bool insert_edge(uintE u, uintE v, W w) {
    Distance ww = edge_weight(w);
    overlay_out[u].push_back({v, ww});

    if (dist[u] == kMaxDistance) return false;

    if (dist[u] > kMaxDistance - ww) return false;

    Distance cand = dist[u] + ww;
    if (cand >= dist[v]) return false;

    auto seeds = sequence<std::pair<uintE, Distance>>(1);
    seeds[0] = std::make_pair(v, cand);

    auto affected = repair_from_seeds(seeds);
    maintain_support_graph_after_repair(affected);

    return !affected.empty();
  }

  bool insert_undirected_edge(uintE u, uintE v, W w) {
    bool c1 = insert_edge(u, v, w);
    bool c2 = insert_edge(v, u, w);
    return c1 || c2;
  }

  bool insert_center(uintE c) {
    if (!virtual_source_attached[c]) {
      virtual_source_attached[c] = true;
      attached_centers.push_back(c);
    }

    if (static_cast<Distance>(0) >= dist[c]) {
      return false;
    }

    auto seeds = sequence<std::pair<uintE, Distance>>(1);
    seeds[0] = std::make_pair(c, static_cast<Distance>(0));

    auto affected = repair_from_seeds(seeds);
    maintain_support_graph_after_repair(affected);

    return !affected.empty();
  }

  inline bool is_support_edge_weighted(uintE u, uintE v, Distance w) const {
    if (dist[u] == kMaxDistance || dist[v] == kMaxDistance) return false;
    if (dist[u] > kMaxDistance - w) return false;
    return dist[u] + w == dist[v];
  }

  inline bool is_support_edge(uintE u, uintE v, W w) const {
    return is_support_edge_weighted(u, v, edge_weight(w));
  }

  std::vector<uintE> compute_support_adj_for_vertex(uintE u) const {
    std::vector<uintE> out;

    auto map_f = [&](const uintE& src, const uintE& v, const W& w) {
      (void)src;
      if (is_support_edge(u, v, w)) {
        out.push_back(v);
      }
    };

    G.get_vertex(u).out_neighbors().map(map_f, false);

    for (const auto& e : overlay_out[u]) {
      if (is_support_edge_weighted(u, e.to, e.w)) {
        out.push_back(e.to);
      }
    }

    return out;
  }

  void rebuild_support_csr_from_adj() {
    const size_t n = G.n;

    support_offsets = sequence<uintE>(n + 1);
    support_offsets[0] = 0;

    for (size_t i = 0; i < n; ++i) {
      support_offsets[i + 1] =
          support_offsets[i] + static_cast<uintE>(support_adj[i].size());
    }

    support_edges = sequence<uintE>(support_offsets[n]);

    parallel_for(0, n, [&](size_t u) {
      size_t idx = support_offsets[u];
      for (uintE v : support_adj[u]) {
        support_edges[idx++] = v;
      }
    });
  }

  void rebuild_support_graph() {
    const size_t n = G.n;

    if (support_adj.size() != n) {
      support_adj = std::vector<std::vector<uintE>>(n);
    }

    parallel_for(0, n, [&](size_t u) {
      support_adj[u] =
          compute_support_adj_for_vertex(static_cast<uintE>(u));
    });

    rebuild_support_csr_from_adj();
  }

  void update_support_graph_around(const sequence<uintE>& affected_vertices) {
    const size_t n = G.n;

    if (affected_vertices.empty()) return;

    auto touched = sequence<uint8_t>(n, static_cast<uint8_t>(0));

    parallel_for(0, affected_vertices.size(), [&](size_t i) {
      uintE u = affected_vertices[i];

      gbbs::atomic_compare_and_swap(&touched[u],
                                    static_cast<uint8_t>(0),
                                    static_cast<uint8_t>(1));

      auto map_f = [&](const uintE& src, const uintE& v, const W& w) {
        (void)src;
        (void)w;

        gbbs::atomic_compare_and_swap(&touched[v],
                                      static_cast<uint8_t>(0),
                                      static_cast<uint8_t>(1));
      };

      G.get_vertex(u).out_neighbors().map(map_f, false);

      for (const auto& e : overlay_out[u]) {
        gbbs::atomic_compare_and_swap(&touched[e.to],
                                      static_cast<uint8_t>(0),
                                      static_cast<uint8_t>(1));
      }
    });

    auto touched_vertices = parlay::filter(
        parlay::tabulate(n, [&](size_t i) {
          return static_cast<uintE>(i);
        }),
        [&](uintE v) { return touched[v] != 0; });

    parallel_for(0, touched_vertices.size(), [&](size_t i) {
      uintE u = touched_vertices[i];
      support_adj[u] = compute_support_adj_for_vertex(u);
    });

    rebuild_support_csr_from_adj();
  }
};

}  // namespace gbbs