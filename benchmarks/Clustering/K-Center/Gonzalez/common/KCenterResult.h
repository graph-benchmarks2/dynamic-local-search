#pragma once

#include <cstddef>

#include "gbbs/gbbs.h"
#include "parlay/sequence.h"

namespace gbbs {

template <class Distance>
struct KCenterResult {
  parlay::sequence<uintE> centers;
  Distance max_dist_to_centers;
  size_t unreachable_vertices;
};

}  // namespace gbbs