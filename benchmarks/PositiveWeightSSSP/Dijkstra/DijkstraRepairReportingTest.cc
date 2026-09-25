#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

#include "benchmarks/PositiveWeightSSSP/Dijkstra/DijkstraRepair.h"
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
  using W = gbbs::uintE;
  using edge = std::tuple<gbbs::uintE, gbbs::uintE, W>;
  using Graph = gbbs::symmetric_graph<gbbs::symmetric_vertex, W>;

  gbbs::sequence<edge> edges = {
      edge{0, 1, 5},
      edge{1, 2, 5},
      edge{2, 3, 5},
      edge{3, 4, 5},
      edge{0, 4, 1},
  };

  Graph G = Graph::from_edges(edges, 5);

  // Distances before inserting weighted shortcut (0,4) of weight 1.
  gbbs::sequence<W> initial = {0, 5, 10, 15, 20};
  gbbs::sequence<std::pair<gbbs::uintE, W>> seeds = {
      {4, 1},
  };

  auto legacy_dist = initial;
  const auto legacy_affected =
      gbbs::DijkstraRepairAffected(G, legacy_dist, seeds);

  auto reported_dist = initial;
  gbbs::DijkstraRepairScratch scratch;
  std::vector<gbbs::uintE> reported_affected;

  gbbs::DijkstraRepair(
    G,
    reported_dist,
    seeds,
    std::numeric_limits<double>::infinity(),
    static_cast<gbbs::AlphaBucketState<gbbs::uintE>*>(nullptr),
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
  gbbs::DijkstraRepair(G, old_api_dist, seeds);
  CheckSequenceEqual(
      old_api_dist,
      legacy_dist,
      "legacy non-reporting DijkstraRepair changed behavior");

  // Reuse scratch and verify that a second repair reports only new changes.
  gbbs::sequence<std::pair<gbbs::uintE, W>> second_seeds = {
      {2, 2},
  };
  reported_affected.clear();

  gbbs::DijkstraRepair(
    G,
    reported_dist,
    second_seeds,
    std::numeric_limits<double>::infinity(),
    static_cast<gbbs::AlphaBucketState<gbbs::uintE>*>(nullptr),
    reported_affected,
    scratch);

  if (reported_affected.size() != 1 || reported_affected[0] != 2) {
    Fail("scratch reuse reported an incorrect affected set");
  }

  std::cout << "DijkstraRepair output-sensitive reporting test passed."
            << std::endl;
  return 0;
}

