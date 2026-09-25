// Sequential Dijkstra implementation in GBBS style.
// Supports positive weighted graphs and also falls back to unit weights
// when the graph is unweighted.

#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "gbbs/gbbs.h"

namespace gbbs {

template <class Graph>
auto Dijkstra(
    Graph& G, uintE src,
    double dist_cap = std::numeric_limits<double>::infinity(),
    bool verbose = true) {
  using W = typename Graph::weight_type;
  using Distance =
      typename std::conditional<std::is_same<W, gbbs::empty>::value, uintE,
                                W>::type;
  constexpr Distance kMaxWeight = std::numeric_limits<Distance>::max();

  const Distance CAP =
      (std::isinf(dist_cap) ? kMaxWeight : static_cast<Distance>(dist_cap));

  const size_t n = G.n;
  auto dists = sequence<Distance>(n, kMaxWeight);

  using PQEntry = std::pair<Distance, uintE>;
  std::priority_queue<PQEntry, std::vector<PQEntry>, std::greater<PQEntry>> pq;

  dists[src] = static_cast<Distance>(0);
  pq.push({static_cast<Distance>(0), src});

  while (!pq.empty()) {
    auto [du, u] = pq.top();
    pq.pop();

    // Skip stale priority-queue entries.
    if (du != dists[u]) continue;
    if (du > CAP) continue;

    auto map_f = [&](const uintE& s, const uintE& d, const W& wgh) {
      Distance w;
      if constexpr (std::is_same<W, gbbs::empty>::value) {
        w = static_cast<Distance>(1);
      } else {
        w = static_cast<Distance>(wgh);
      }

      if (dists[s] > kMaxWeight - w) return;

      Distance nd = dists[s] + w;
      if (nd > CAP) return;

      if (nd < dists[d]) {
        dists[d] = nd;
        pq.push({nd, d});
      }
    };

    G.get_vertex(u).out_neighbors().map(map_f, false);
  }

  auto get_dist = [&](size_t i) -> Distance {
    return (dists[i] == kMaxWeight) ? static_cast<Distance>(0) : dists[i];
  };
  auto dist_im = parlay::delayed_seq<Distance>(n, get_dist);
  if (verbose) {
    std::cout << "max_dist = " << parlay::reduce_max(dist_im) << std::endl;
  };

  return dists;
}

}  // namespace gbbs