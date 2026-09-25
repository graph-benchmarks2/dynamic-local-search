// This file provides a parallel BFS variant for unweighted graphs whose primary goal is to
// maintain and return an array of distances, not a full BFS tree.

#pragma once

#include <cmath>
#include <limits>

#include "gbbs/gbbs.h"

namespace gbbs {

template <class W>
struct BFSDist_F {
  uintE* dists;
  uintE next_dist;
  uintE cap;

  BFSDist_F(uintE* _dists, uintE _next_dist, uintE _cap)
      : dists(_dists), next_dist(_next_dist), cap(_cap) {}

  inline bool update(uintE s, uintE d, W w) {
    (void)s;
    (void)w;
    if (next_dist > cap) return false;
    if (dists[d] == UINT_E_MAX) {
      dists[d] = next_dist;
      return true;
    }
    return false;
  }

  inline bool updateAtomic(uintE s, uintE d, W w) {
    (void)s;
    (void)w;
    if (next_dist > cap) return false;
    return gbbs::atomic_compare_and_swap(&dists[d], UINT_E_MAX, next_dist);
  }

  inline bool cond(uintE d) const {
    return (next_dist <= cap && dists[d] == UINT_E_MAX);
  }
};

template <class Graph>
inline sequence<uintE> BFSdist(
    Graph& G, uintE src,
    double dist_cap = std::numeric_limits<double>::infinity(),
    bool verbose = true) {
  using W = typename Graph::weight_type;

  const uintE CAP =
      (std::isinf(dist_cap) ? UINT_E_MAX : static_cast<uintE>(dist_cap));

  auto dists =
      sequence<uintE>::from_function(G.n, [&](size_t i) { return UINT_E_MAX; });
  dists[src] = 0;

  vertexSubset Frontier(G.n, src);
  uintE cur_dist = 0;

  while (!Frontier.isEmpty()) {
    if (cur_dist >= CAP) break;

    ++cur_dist;
    Frontier = edgeMap(
        G, Frontier, BFSDist_F<W>(dists.begin(), cur_dist, CAP), -1,
        sparse_blocked | dense_parallel);
  }

  auto dist_im = parlay::delayed_seq<uintE>(
      G.n, [&](size_t i) { return (dists[i] == UINT_E_MAX) ? 0 : dists[i]; });
  if (verbose) {
    std::cout << "max_dist = " << parlay::reduce_max(dist_im) << std::endl;
  }

  return dists;
}

}  // namespace gbbs