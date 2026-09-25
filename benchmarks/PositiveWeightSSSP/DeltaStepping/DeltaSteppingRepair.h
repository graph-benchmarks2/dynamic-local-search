// This file provides a parallel implementation of a Delta Stepping variant
// that can use a previously computed distance array.
//
// Optionally, multiplicative alpha-bucket bookkeeping can be enabled
// by passing a pointer. In that mode, we snapshot the
// old distances, mark changed vertices during the distance update run and update
// bucket memberships after.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "gbbs/gbbs.h"
#include "benchmarks/PositiveWeightSSSP/common/AlphaBucketState.h"

namespace gbbs {

// Repair vertex

template <class W, class Distance>
struct Repair_Visit_F {
  sequence<std::pair<Distance, bool>>& dists;
  Distance cap;
  sequence<uint8_t>* changed;

  Repair_Visit_F(sequence<std::pair<Distance, bool>>& _dists, Distance _cap,
                 sequence<uint8_t>* _changed = nullptr)
      : dists(_dists), cap(_cap), changed(_changed) {}

  inline std::optional<Distance> update(const uintE& s, const uintE& d,
                                        const W& w) {
    Distance cur = dists[d].first;
    Distance cand;
    if constexpr (std::is_same<W, gbbs::empty>::value) {
      cand = dists[s].first + 1;
    } else {
      cand = dists[s].first + w;
    }

    if (cand > cap) return std::nullopt;

    if (cand < cur) {
      if (changed != nullptr) {
        (*changed)[d] = 1;
      }
      if (!dists[d].second) {
        dists[d] = {cand, true};
        return std::optional<Distance>(cur);
      }
      dists[d].first = cand;
    }
    return std::nullopt;
  }

  inline std::optional<Distance> updateAtomic(const uintE& s, const uintE& d,
                                              const W& w) {
    Distance cur = dists[d].first;
    Distance cand;
    if constexpr (std::is_same<W, gbbs::empty>::value) {
      cand = dists[s].first + 1;
    } else {
      cand = dists[s].first + w;
    }

    if (cand > cap) return std::nullopt;

    if (cand < cur) {
      gbbs::write_min(&(dists[d].first), cand);

      if (changed != nullptr) {
        gbbs::atomic_compare_and_swap(&((*changed)[d]), static_cast<uint8_t>(0),
                                      static_cast<uint8_t>(1));
      }

      if (!dists[d].second &&
          gbbs::atomic_compare_and_swap(&dists[d].second, false, true)) {
        return std::optional<Distance>(cur);
      }
    }
    return std::nullopt;
  }

  inline bool cond(const uintE& d) const { return true; }
};

// Internal implementation of DeltaSteppingRepair that does not use buckets

template <class Graph, class Distance>
void DeltaSteppingRepairImpl(
    Graph& G, sequence<Distance>& dist,
    const sequence<std::pair<uintE, Distance>>& seeds, double delta,
    size_t num_buckets = 128,
    double dist_cap = std::numeric_limits<double>::infinity(),
    sequence<uint8_t>* changed = nullptr) {
  using W = typename Graph::weight_type;
  constexpr Distance kMaxWeight = std::numeric_limits<Distance>::max();
  (void)delta;
  (void)num_buckets;

  const Distance CAP =
      (std::isinf(dist_cap) ? kMaxWeight : static_cast<Distance>(dist_cap));

  const size_t n = G.n;

  if (dist.size() != n) {
    std::cout << "DeltaSteppingRepair: dist.size() != G.n" << std::endl;
    exit(-1);
  }

  if (changed != nullptr && changed->size() != n) {
    std::cout << "DeltaSteppingRepair: changed.size() != G.n" << std::endl;
    exit(-1);
  }

  auto dists = sequence<std::pair<Distance, bool>>::from_function(
      n, [&](size_t i) { return std::make_pair(dist[i], false); });

  parallel_for(0, seeds.size(), [&](size_t i) {
    const auto& [v, new_dist] = seeds[i];
    if (new_dist > CAP) return;
    if (new_dist < dists[v].first) {
      gbbs::write_min(&(dists[v].first), new_dist);
      gbbs::atomic_compare_and_swap(&dists[v].second, false, true);
      if (changed != nullptr) {
        gbbs::atomic_compare_and_swap(&((*changed)[v]), static_cast<uint8_t>(0),
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

    auto res = edgeMapData<Distance>(
        G, active_vs, Repair_Visit_F<W, Distance>(dists, CAP, changed), G.m / 20,
        fl);

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

template <class Graph, class Distance>
sequence<uintE> DeltaSteppingRepairAffected(
    Graph& G,
    sequence<Distance>& dist,
    const sequence<std::pair<uintE, Distance>>& seeds,
    double delta,
    size_t num_buckets = 128,
    double dist_cap = std::numeric_limits<double>::infinity()) {
  const size_t n = G.n;
  auto changed = sequence<uint8_t>(n, static_cast<uint8_t>(0));

  DeltaSteppingRepairImpl(G, dist, seeds, delta, num_buckets, dist_cap,
                          &changed);

  return parlay::filter(
      parlay::tabulate(n, [&](size_t i) { return static_cast<uintE>(i); }),
      [&](uintE v) { return changed[v] != 0; });
}

template <class Graph, class Distance>
void DeltaSteppingRepair(
    Graph& G, sequence<Distance>& dist,
    const sequence<std::pair<uintE, Distance>>& seeds, double delta,
    size_t num_buckets = 128,
    double dist_cap = std::numeric_limits<double>::infinity(),
    AlphaBucketState<Distance>* alpha_state = nullptr) {
  const size_t n = G.n;

  if (alpha_state == nullptr) {
    DeltaSteppingRepairImpl(G, dist, seeds, delta, num_buckets, dist_cap,
                            nullptr);
    return;
  }

  auto old_dist = sequence<Distance>::from_function(n, [&](size_t i) {
    return dist[i];
  });
  auto changed = sequence<uint8_t>(n, static_cast<uint8_t>(0));

  DeltaSteppingRepairImpl(G, dist, seeds, delta, num_buckets, dist_cap,
                          &changed);

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

template <class Graph, class Distance>
void DeltaSteppingRepair(
    Graph& G, sequence<Distance>& dist,
    const sequence<uintE>& active_vertices, double delta,
    size_t num_buckets = 128,
    double dist_cap = std::numeric_limits<double>::infinity(),
    AlphaBucketState<Distance>* alpha_state = nullptr) {
  constexpr Distance kMaxWeight = std::numeric_limits<Distance>::max();

  auto seeds = sequence<std::pair<uintE, Distance>>::from_function(
      active_vertices.size(), [&](size_t i) {
        uintE v = active_vertices[i];
        return std::make_pair(v, dist[v]);
      });

  parallel_for(0, active_vertices.size(), [&](size_t i) {
    uintE v = active_vertices[i];
    dist[v] = kMaxWeight;
  });

  DeltaSteppingRepair(G, dist, seeds, delta, num_buckets, dist_cap, alpha_state);
}

template <class Graph, class Distance>
void DeltaSteppingRepair(
    Graph& G, sequence<Distance>& dist, uintE src, double delta,
    size_t num_buckets = 128,
    double dist_cap = std::numeric_limits<double>::infinity(),
    AlphaBucketState<Distance>* alpha_state = nullptr) {
  auto active = sequence<uintE>(1, src);
  DeltaSteppingRepair(G, dist, active, delta, num_buckets, dist_cap,
                      alpha_state);
}

}  // namespace gbbs