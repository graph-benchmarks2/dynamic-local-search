#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

#include "benchmarks/BFS/BFSdist/BFSdistRepairSeq.h"
#include "gbbs/gbbs.h"

namespace {

void Fail(const char* message) {
  std::cout << "ERROR: " << message << std::endl;
  std::exit(-1);
}

template <class T>
void CheckSequenceEqual(
    const gbbs::sequence<T>& a,
    const gbbs::sequence<T>& b,
    const char* message) {
  if (a.size() != b.size()) Fail(message);
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i] != b[i]) Fail(message);
  }
}

std::vector<gbbs::uintE> Sorted(std::vector<gbbs::uintE> values) {
  std::sort(values.begin(), values.end());
  return values;
}

std::vector<gbbs::uintE> Sorted(const gbbs::sequence<gbbs::uintE>& values) {
  std::vector<gbbs::uintE> result(values.begin(), values.end());
  std::sort(result.begin(), result.end());
  return result;
}

}  // namespace

int main() {
  using W = gbbs::empty;
  using edge = std::tuple<gbbs::uintE, gbbs::uintE, W>;
  using Graph = gbbs::symmetric_graph<gbbs::symmetric_vertex, W>;

  gbbs::sequence<edge> edges = {
      edge{0, 1, W{}},
      edge{1, 2, W{}},
      edge{2, 3, W{}},
      edge{3, 4, W{}},
      edge{0, 4, W{}},
  };

  Graph G = Graph::from_edges(edges, 5);

  // Distances before inserting shortcut (0,4).
  gbbs::sequence<gbbs::uintE> initial = {0, 1, 2, 3, 4};
  gbbs::sequence<std::pair<gbbs::uintE, gbbs::uintE>> seeds = {
      {4, 1},
  };

  auto legacy_dist = initial;
  const auto legacy_affected =
      gbbs::BFSdistRepairSeqAffected(G, legacy_dist, seeds);

  auto reported_dist = initial;
  gbbs::BFSdistRepairSeqScratch scratch;
  std::vector<gbbs::uintE> reported_affected;

  gbbs::BFSdistRepairSeq(
      G,
      reported_dist,
      seeds,
      std::numeric_limits<double>::infinity(),
      nullptr,
      reported_affected,
      scratch);

  CheckSequenceEqual(
      reported_dist,
      legacy_dist,
      "reported repair distances disagree with legacy Affected repair");

  if (Sorted(reported_affected) != Sorted(legacy_affected)) {
    Fail("reported affected vertices disagree with legacy Affected repair");
  }

  // The old non-reporting call path must remain semantically unchanged.
  auto old_api_dist = initial;
  gbbs::BFSdistRepairSeq(G, old_api_dist, seeds);
  CheckSequenceEqual(
      old_api_dist,
      legacy_dist,
      "legacy non-reporting BFSdistRepairSeq changed behavior");

  // Reuse the same scratch object to make sure epochs do not leak state.
  gbbs::sequence<std::pair<gbbs::uintE, gbbs::uintE>> second_seeds = {
      {2, 1},
  };
  reported_affected.clear();

  gbbs::BFSdistRepairSeq(
      G,
      reported_dist,
      second_seeds,
      std::numeric_limits<double>::infinity(),
      nullptr,
      reported_affected,
      scratch);

  if (reported_affected.size() != 1 || reported_affected[0] != 2) {
    Fail("scratch reuse reported an incorrect affected set");
  }

  std::cout << "BFSdistRepairSeq output-sensitive reporting test passed."
            << std::endl;
  return 0;
}

