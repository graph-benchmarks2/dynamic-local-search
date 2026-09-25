#pragma once

#include "gbbs/gbbs.h"

namespace gbbs {
namespace kmedian {

struct DistanceDecrease {
    uintE p = 0;
    uintE q = 0;
    double new_distance = 0.0;
};

}  // namespace kmedian
}  // namespace gbbs
