// This file provides a sequential BFS variant for unweighted graphs whose
// primary goal is to maintain and return an array of distances, not a full BFS
// tree. It is intended for deliberate single-core usage.

#pragma once

#include <cmath>
#include <limits>
#include <queue>

#include "gbbs/gbbs.h"

namespace gbbs {

template <class Graph>
inline sequence<uintE> BFSdistSeq(
    Graph& G, uintE src,
    double dist_cap = std::numeric_limits<double>::infinity(),
    bool verbose = true) {
  const uintE CAP =
      (std::isinf(dist_cap) ? UINT_E_MAX : static_cast<uintE>(dist_cap));

  auto dist =
      sequence<uintE>::from_function(G.n, [&](size_t i) { return UINT_E_MAX; });

  if (G.n == 0) return dist;

  dist[src] = 0;

  std::queue<uintE> q;
  q.push(src);

  while (!q.empty()) {
    uintE u = q.front();
    q.pop();

    if (dist[u] >= CAP) continue;

    uintE cand = dist[u] + 1;

    auto map_f = [&](const uintE& src, const uintE& v,
                     const typename Graph::weight_type& w) {
      (void)src;
      (void)w;

      if (cand <= CAP && cand < dist[v]) {
        dist[v] = cand;
        q.push(v);
      }
    };

    G.get_vertex(u).out_neighbors().map(map_f, false);
  }

  uintE max_dist = 0;
  for (size_t i = 0; i < G.n; ++i) {
    if (dist[i] != UINT_E_MAX) {
      max_dist = std::max(max_dist, dist[i]);
    }
  }

  if (verbose) {
    std::cout << "max_dist = " << max_dist << std::endl;
  }

  return dist;
}

}  // namespace gbbs
