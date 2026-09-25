// Shared state for distance-based bucketing.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "gbbs/gbbs.h"

namespace gbbs {

template <class Distance>
struct AlphaBucketState {
  static constexpr size_t kNoBucket = std::numeric_limits<size_t>::max();

  size_t n;
  double alpha;
  size_t B;
  Distance inf_dist;

  std::vector<std::vector<uintE>> H;

  sequence<size_t> bucket_of;
  sequence<size_t> pos_in_bucket;
  sequence<uint8_t> is_center;

  AlphaBucketState(size_t _n, double _alpha, size_t _B, Distance _inf_dist)
      : n(_n),
        alpha(_alpha),
        B(_B),
        inf_dist(_inf_dist),
        H(_B + 1),
        bucket_of(sequence<size_t>(_n, kNoBucket)),
        pos_in_bucket(sequence<size_t>(_n, 0)),
        is_center(sequence<uint8_t>(_n, 0)) {}

  inline size_t distance_to_bucket(Distance d) const {
    if (d == static_cast<Distance>(0)) return kNoBucket;
    if (d == inf_dist) return B;

    double x = static_cast<double>(d);
    double idx = std::ceil(std::log(x) / std::log(alpha));

    if (idx < 1.0) return 1;
    size_t b = static_cast<size_t>(idx);
    return std::min(b, B);
  }

  inline void insert_vertex(uintE v, size_t b) {
    if (b == kNoBucket) return;
    H[b].push_back(v);
    bucket_of[v] = b;
    pos_in_bucket[v] = H[b].size() - 1;
  }

  inline void remove_vertex(uintE v) {
    size_t b = bucket_of[v];
    if (b == kNoBucket) return;

    auto& bucket = H[b];
    size_t pos = pos_in_bucket[v];
    uintE last = bucket.back();

    bucket[pos] = last;
    pos_in_bucket[last] = pos;
    bucket.pop_back();

    bucket_of[v] = kNoBucket;
  }

  inline void move_vertex(uintE v, size_t new_bucket) {
    size_t old_bucket = bucket_of[v];
    if (old_bucket == new_bucket) return;
    remove_vertex(v);
    insert_vertex(v, new_bucket);
  }

  inline void update_vertex_distance(uintE v, Distance new_dist) {
    if (is_center[v]) return;
    size_t new_bucket = distance_to_bucket(new_dist);
    move_vertex(v, new_bucket);
  }

  inline void set_center(uintE v) {
    if (is_center[v]) return;
    remove_vertex(v);
    is_center[v] = 1;
  }

  inline bool bucket_empty(size_t b) const { return H[b].empty(); }

  inline void initialize_all_unreachable() {
    for (size_t b = 1; b <= B; ++b) H[b].clear();
    for (size_t i = 0; i < n; ++i) {
      is_center[i] = 0;
      bucket_of[i] = kNoBucket;
      pos_in_bucket[i] = 0;
    }
    H[B].reserve(n);
    for (uintE v = 0; v < n; ++v) {
      insert_vertex(v, B);
    }
  }
};

}  // namespace gbbs