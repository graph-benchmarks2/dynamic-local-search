// This file provides a sequential BFS repair variant for unweighted graphs whose
// primary goal is to maintain and return an array of distances, not a full BFS
// tree. It is intended for deliberate single-core usage.
// In addition, it can use a previously computed distance array.
//
// Optionally, multiplicative alpha-bucket bookkeeping can be enabled by passing
// a pointer. In that mode, we snapshot the old distances, mark changed vertices
// during the distance update run and update bucket memberships after.

#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "gbbs/gbbs.h"
#include "benchmarks/PositiveWeightSSSP/common/AlphaBucketState.h"

namespace gbbs {

// Reusable scratch space for the optional output-sensitive reporting overload.
// Existing BFSdistRepairSeq call sites do not allocate or use this state.
struct BFSdistRepairSeqScratch {
  std::vector<uint32_t> in_queue_epoch;
  std::vector<uint32_t> reported_epoch;
  uint32_t epoch = 0;

  uint32_t begin(size_t n) {
    if (in_queue_epoch.size() != n) {
      in_queue_epoch.assign(n, 0);
      reported_epoch.assign(n, 0);
      epoch = 1;
      return epoch;
    }

    ++epoch;
    if (epoch == 0) {
      std::fill(in_queue_epoch.begin(), in_queue_epoch.end(), 0);
      std::fill(reported_epoch.begin(), reported_epoch.end(), 0);
      epoch = 1;
    }
    return epoch;
  }
};

template <class Graph>
void BFSdistRepairSeqImpl(
    Graph& G, sequence<uintE>& dist,
    const sequence<std::pair<uintE, uintE>>& seeds,
    double dist_cap = std::numeric_limits<double>::infinity(),
    sequence<uint8_t>* changed = nullptr) {
  const uintE CAP =
      (std::isinf(dist_cap) ? UINT_E_MAX : static_cast<uintE>(dist_cap));

  const size_t n = G.n;

  if (dist.size() != n) {
    std::cout << "BFSdistRepairSeq: dist.size() != G.n" << std::endl;
    exit(-1);
  }

  if (changed != nullptr && changed->size() != n) {
    std::cout << "BFSdistRepairSeq: changed.size() != G.n" << std::endl;
    exit(-1);
  }

  std::queue<uintE> q;
  std::vector<uint8_t> in_queue(n, static_cast<uint8_t>(0));

  for (size_t i = 0; i < seeds.size(); ++i) {
    const auto& [v, new_dist] = seeds[i];

    if (v >= n) continue;
    if (new_dist > CAP) continue;

    if (new_dist < dist[v]) {
      dist[v] = new_dist;

      if (changed != nullptr) {
        (*changed)[v] = 1;
      }

      if (!in_queue[v]) {
        q.push(v);
        in_queue[v] = 1;
      }
    }
  }

  while (!q.empty()) {
    uintE u = q.front();
    q.pop();
    in_queue[u] = 0;

    if (dist[u] >= CAP) continue;

    uintE cand = dist[u] + 1;

    auto map_f = [&](const uintE& src, const uintE& v,
                     const typename Graph::weight_type& w) {
      (void)src;
      (void)w;

      if (cand <= CAP && cand < dist[v]) {
        dist[v] = cand;

        if (changed != nullptr) {
          (*changed)[v] = 1;
        }

        if (!in_queue[v]) {
          q.push(v);
          in_queue[v] = 1;
        }
      }
    };

    G.get_vertex(u).out_neighbors().map(map_f, false);
  }
}

template <class Graph>
sequence<uintE> BFSdistRepairSeqAffected(
    Graph& G,
    sequence<uintE>& dist,
    const sequence<std::pair<uintE, uintE>>& seeds,
    double dist_cap = std::numeric_limits<double>::infinity()) {
  const size_t n = G.n;
  auto changed = sequence<uint8_t>(n, static_cast<uint8_t>(0));

  BFSdistRepairSeqImpl(G, dist, seeds, dist_cap, &changed);

  return parlay::filter(
      parlay::tabulate(n, [&](size_t i) { return static_cast<uintE>(i); }),
      [&](uintE v) { return changed[v] != 0; });
}

template <class Graph>
void BFSdistRepairSeq(
    Graph& G, sequence<uintE>& dist,
    const sequence<std::pair<uintE, uintE>>& seeds,
    double dist_cap = std::numeric_limits<double>::infinity(),
    AlphaBucketState<uintE>* alpha_state = nullptr) {
  const size_t n = G.n;

  if (alpha_state == nullptr) {
    BFSdistRepairSeqImpl(G, dist, seeds, dist_cap, nullptr);
    return;
  }

  auto old_dist = sequence<uintE>::from_function(n, [&](size_t i) {
    return dist[i];
  });
  auto changed = sequence<uint8_t>(n, static_cast<uint8_t>(0));

  BFSdistRepairSeqImpl(G, dist, seeds, dist_cap, &changed);

  for (uintE v = 0; v < n; ++v) {
    if (changed[v] == 0) continue;
    if (alpha_state->is_center[v]) continue;
    if (dist[v] < old_dist[v]) {
      alpha_state->update_vertex_distance(v, dist[v]);
    }
  }
}

// Output-sensitive reporting overload. It preserves the existing repair
// semantics but directly appends each vertex whose distance decreases, exactly
// once, to affected_vertices. Reusable scratch avoids O(n) initialization and
// O(n) filtering on every repair. Existing overloads above are unchanged.
template <class Graph>
void BFSdistRepairSeq(
    Graph& G, sequence<uintE>& dist,
    const sequence<std::pair<uintE, uintE>>& seeds,
    double dist_cap,
    AlphaBucketState<uintE>* alpha_state,
    std::vector<uintE>& affected_vertices,
    BFSdistRepairSeqScratch& scratch) {
  const uintE CAP =
      (std::isinf(dist_cap) ? UINT_E_MAX : static_cast<uintE>(dist_cap));

  const size_t n = G.n;
  if (dist.size() != n) {
    std::cout << "BFSdistRepairSeq: dist.size() != G.n" << std::endl;
    exit(-1);
  }

  affected_vertices.clear();
  const uint32_t epoch = scratch.begin(n);
  std::queue<uintE> q;

  auto report = [&](uintE v) {
    if (scratch.reported_epoch[v] != epoch) {
      scratch.reported_epoch[v] = epoch;
      affected_vertices.push_back(v);
    }
  };

  auto enqueue = [&](uintE v) {
    if (scratch.in_queue_epoch[v] != epoch) {
      scratch.in_queue_epoch[v] = epoch;
      q.push(v);
    }
  };

  for (size_t i = 0; i < seeds.size(); ++i) {
    const auto& [v, new_dist] = seeds[i];
    if (v >= n || new_dist > CAP) continue;

    if (new_dist < dist[v]) {
      dist[v] = new_dist;
      report(v);
      enqueue(v);
    }
  }

  while (!q.empty()) {
    const uintE u = q.front();
    q.pop();
    scratch.in_queue_epoch[u] = 0;

    if (dist[u] >= CAP) continue;
    const uintE cand = dist[u] + 1;

    auto map_f = [&](const uintE& src, const uintE& v,
                     const typename Graph::weight_type& w) {
      (void)src;
      (void)w;

      if (cand <= CAP && cand < dist[v]) {
        dist[v] = cand;
        report(v);
        enqueue(v);
      }
    };

    G.get_vertex(u).out_neighbors().map(map_f, false);
  }

  if (alpha_state != nullptr) {
    for (uintE v : affected_vertices) {
      if (!alpha_state->is_center[v]) {
        alpha_state->update_vertex_distance(v, dist[v]);
      }
    }
  }
}

template <class Graph>
void BFSdistRepairSeq(
    Graph& G, sequence<uintE>& dist,
    const sequence<uintE>& active_vertices,
    double dist_cap = std::numeric_limits<double>::infinity(),
    AlphaBucketState<uintE>* alpha_state = nullptr) {
  auto seeds = sequence<std::pair<uintE, uintE>>::from_function(
      active_vertices.size(), [&](size_t i) {
        uintE v = active_vertices[i];
        return std::make_pair(v, dist[v]);
      });

  for (size_t i = 0; i < active_vertices.size(); ++i) {
    uintE v = active_vertices[i];
    dist[v] = UINT_E_MAX;
  }

  BFSdistRepairSeq(G, dist, seeds, dist_cap, alpha_state);
}

template <class Graph>
void BFSdistRepairSeq(
    Graph& G, sequence<uintE>& dist, uintE src,
    double dist_cap = std::numeric_limits<double>::infinity(),
    AlphaBucketState<uintE>* alpha_state = nullptr) {
  auto active = sequence<uintE>(1, src);
  BFSdistRepairSeq(G, dist, active, dist_cap, alpha_state);
}

}  // namespace gbbs
