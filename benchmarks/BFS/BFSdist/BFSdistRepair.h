// This file provides a parallel implementation of a BFS variant for unweighted graphs whose
// primary goal is to maintain and return an array of distances, not a full BFS tree.
// In addition, it can use a previously computed distance array.
//
// Optionally, multiplicative alpha-bucket bookkeeping can be enabled
// by passing a pointer. In that mode, we snapshot the
// old distances, mark changed vertices during the distance update run and update
// bucket memberships after.

#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <tuple>
#include <utility>

#include "gbbs/gbbs.h"
#include "gbbs/edge_map_data.h"
#include "gbbs/vertex_subset.h"
#include "benchmarks/PositiveWeightSSSP/common/AlphaBucketState.h"

namespace gbbs {

template <class W>
struct BFSDistRepair_F {
  sequence<std::pair<uintE, bool>>& dists;
  uintE cap;
  sequence<uint8_t>* changed;

  BFSDistRepair_F(sequence<std::pair<uintE, bool>>& _dists, uintE _cap,
                  sequence<uint8_t>* _changed = nullptr)
      : dists(_dists), cap(_cap), changed(_changed) {}

  inline std::optional<uintE> update(const uintE& s, const uintE& d, const W& w) {
    (void)w;
    uintE cur = dists[d].first;
    uintE cand = dists[s].first + 1;

    if (cand > cap) return std::nullopt;

    if (cand < cur) {
      if (changed != nullptr) {
        (*changed)[d] = 1;
      }
      if (!dists[d].second) {
        dists[d] = {cand, true};
        return std::optional<uintE>(cur);
      }
      dists[d].first = cand;
    }
    return std::nullopt;
  }

  inline std::optional<uintE> updateAtomic(const uintE& s, const uintE& d,
                                           const W& w) {
    (void)w;
    uintE cur = dists[d].first;
    uintE cand = dists[s].first + 1;

    if (cand > cap) return std::nullopt;

    if (cand < cur) {
      gbbs::write_min(&(dists[d].first), cand);

      if (changed != nullptr) {
        gbbs::atomic_compare_and_swap(&((*changed)[d]),
                                      static_cast<uint8_t>(0),
                                      static_cast<uint8_t>(1));
      }

      if (!dists[d].second &&
          gbbs::atomic_compare_and_swap(&dists[d].second, false, true)) {
        return std::optional<uintE>(cur);
      }
    }
    return std::nullopt;
  }

  inline bool cond(const uintE& d) const { return true; }
};

template <class Graph>
void BFSdistRepairImpl(
    Graph& G, sequence<uintE>& dist,
    const sequence<std::pair<uintE, uintE>>& seeds,
    double dist_cap = std::numeric_limits<double>::infinity(),
    sequence<uint8_t>* changed = nullptr) {
  using W = typename Graph::weight_type;

  const uintE CAP =
      (std::isinf(dist_cap) ? UINT_E_MAX : static_cast<uintE>(dist_cap));

  const size_t n = G.n;

  if (dist.size() != n) {
    std::cout << "BFSdistRepair: dist.size() != G.n" << std::endl;
    exit(-1);
  }

  if (changed != nullptr && changed->size() != n) {
    std::cout << "BFSdistRepair: changed.size() != G.n" << std::endl;
    exit(-1);
  }

  auto dists = sequence<std::pair<uintE, bool>>::from_function(
      n, [&](size_t i) { return std::make_pair(dist[i], false); });

  // Apply seeded decreases.
  parallel_for(0, seeds.size(), [&](size_t i) {
    const auto& [v, new_dist] = seeds[i];
    if (new_dist > CAP) return;
    if (new_dist < dists[v].first) {
      gbbs::write_min(&(dists[v].first), new_dist);
      gbbs::atomic_compare_and_swap(&dists[v].second, false, true);
      if (changed != nullptr) {
        gbbs::atomic_compare_and_swap(&((*changed)[v]),
                                      static_cast<uint8_t>(0),
                                      static_cast<uint8_t>(1));
      }
    }
  });

  sequence<uintE> active = parlay::filter(
      parlay::tabulate(n, [&](size_t i) { return static_cast<uintE>(i); }),
      [&](uintE v) { return dists[v].second; });

  flags fl = no_dense;

  while (active.size() > 0) {
    parallel_for(0, active.size(), [&](size_t i) {
      uintE v = active[i];
      dists[v].second = false;
    });

    auto active_vs = vertexSubset(n, std::move(active));

    auto res = edgeMapData<uintE>(
        G, active_vs, BFSDistRepair_F<W>(dists, CAP, changed), G.m / 20, fl);

    if (res.dense()) {
      active = parlay::filter(
          parlay::tabulate(n, [&](size_t i) { return static_cast<uintE>(i); }),
          [&](uintE v) { return std::get<0>(res.d[v]); });
    } else {
      active = sequence<uintE>::from_function(
          res.size(), [&](size_t i) { return std::get<0>(res.s[i]); });
    }
  }

  parallel_for(0, n, [&](size_t i) { dist[i] = dists[i].first; });
}

template <class Graph>
sequence<uintE> BFSdistRepairAffected(
    Graph& G,
    sequence<uintE>& dist,
    const sequence<std::pair<uintE, uintE>>& seeds,
    double dist_cap = std::numeric_limits<double>::infinity()) {
  const size_t n = G.n;
  auto changed = sequence<uint8_t>(n, static_cast<uint8_t>(0));

  BFSdistRepairImpl(G, dist, seeds, dist_cap, &changed);

  return parlay::filter(
      parlay::tabulate(n, [&](size_t i) { return static_cast<uintE>(i); }),
      [&](uintE v) { return changed[v] != 0; });
}

template <class Graph>
void BFSdistRepair(
    Graph& G, sequence<uintE>& dist,
    const sequence<std::pair<uintE, uintE>>& seeds,
    double dist_cap = std::numeric_limits<double>::infinity(),
    AlphaBucketState<uintE>* alpha_state = nullptr) {
  const size_t n = G.n;

  if (alpha_state == nullptr) {
    BFSdistRepairImpl(G, dist, seeds, dist_cap, nullptr);
    return;
  }

  auto old_dist = sequence<uintE>::from_function(n, [&](size_t i) {
    return dist[i];
  });
  auto changed = sequence<uint8_t>(n, static_cast<uint8_t>(0));

  BFSdistRepairImpl(G, dist, seeds, dist_cap, &changed);

  auto changed_vertices = parlay::filter(
      parlay::tabulate(n, [&](size_t i) { return static_cast<uintE>(i); }),
      [&](uintE v) { return changed[v] != 0; });

  for (size_t i = 0; i < changed_vertices.size(); ++i) {
    uintE v = changed_vertices[i];
    if (alpha_state->is_center[v]) continue;
    if (dist[v] < old_dist[v]) {
      alpha_state->update_vertex_distance(v, dist[v]);
    }
  }
}

template <class Graph>
void BFSdistRepair(
    Graph& G, sequence<uintE>& dist,
    const sequence<uintE>& active_vertices,
    double dist_cap = std::numeric_limits<double>::infinity(),
    AlphaBucketState<uintE>* alpha_state = nullptr) {
  auto seeds = sequence<std::pair<uintE, uintE>>::from_function(
      active_vertices.size(), [&](size_t i) {
        uintE v = active_vertices[i];
        return std::make_pair(v, dist[v]);
      });

  parallel_for(0, active_vertices.size(), [&](size_t i) {
    uintE v = active_vertices[i];
    dist[v] = UINT_E_MAX;
  });

  BFSdistRepair(G, dist, seeds, dist_cap, alpha_state);
}

template <class Graph>
void BFSdistRepair(
    Graph& G, sequence<uintE>& dist, uintE src,
    double dist_cap = std::numeric_limits<double>::infinity(),
    AlphaBucketState<uintE>* alpha_state = nullptr) {
  auto active = sequence<uintE>(1, src);
  BFSdistRepair(G, dist, active, dist_cap, alpha_state);
}

}  // namespace gbbs