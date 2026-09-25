#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include "gbbs/gbbs.h"

#include "benchmarks/BFS/BFSdist/BFSdistSeq.h"
#include "benchmarks/Clustering/K-Median/common/IncrementalAPSP.h"
#include "benchmarks/PositiveWeightSSSP/Dijkstra/Dijkstra.h"

namespace {

template <class T>
void CheckEqual(
    const T& got,
    const T& expected,
    const char* message) {
  if (!(got == expected)) {
    std::cout << "ERROR: " << message
              << ": expected " << expected
              << ", got " << got
              << std::endl;
    std::exit(-1);
  }
}

template <class Graph>
using APSP = gbbs::kmedian::IncrementalAPSP<Graph>;

template <class Graph>
std::vector<std::vector<typename APSP<Graph>::Distance>>
FreshAPSP(Graph& G) {
  using W = typename Graph::weight_type;
  using Distance = typename APSP<Graph>::Distance;

  std::vector<std::vector<Distance>> result(
      G.n,
      std::vector<Distance>(G.n));

  for (gbbs::uintE source = 0; source < G.n; ++source) {
    gbbs::sequence<Distance> row;

    if constexpr (std::is_same<W, gbbs::empty>::value) {
      row = gbbs::BFSdistSeq(
          G, source, std::numeric_limits<double>::infinity(), /*verbose=*/false);
    } else {
      row = gbbs::Dijkstra(
          G, source, std::numeric_limits<double>::infinity(), /*verbose=*/false);
    }

    for (size_t target = 0; target < G.n; ++target) {
      result[source][target] = row[target];
    }
  }

  return result;
}

template <class Graph>
void VerifyMatrix(
    const APSP<Graph>& apsp,
    const std::vector<std::vector<typename APSP<Graph>::Distance>>& expected,
    const char* message) {
  for (gbbs::uintE u = 0; u < apsp.n(); ++u) {
    for (gbbs::uintE v = 0; v < apsp.n(); ++v) {
      if (apsp.distance(u, v) != expected[u][v]) {
        std::cout << "ERROR: " << message
                  << " at (" << u << "," << v << ")"
                  << ": expected " << expected[u][v]
                  << ", got " << apsp.distance(u, v)
                  << std::endl;
        std::exit(-1);
      }
    }
  }
}

template <class Distance>
std::vector<std::tuple<gbbs::uintE, gbbs::uintE, Distance>>
ExpectedDecreases(
    const std::vector<std::vector<Distance>>& old_dist,
    const std::vector<std::vector<Distance>>& new_dist) {
  std::vector<std::tuple<gbbs::uintE, gbbs::uintE, Distance>> out;

  for (gbbs::uintE u = 0; u < old_dist.size(); ++u) {
    for (gbbs::uintE v = u + 1; v < old_dist.size(); ++v) {
      if (new_dist[u][v] < old_dist[u][v]) {
        out.emplace_back(u, v, new_dist[u][v]);
      }
    }
  }

  return out;
}

template <class Graph>
void VerifyBatch(
    const gbbs::sequence<gbbs::kmedian::DistanceDecrease>& batch,
    const std::vector<std::tuple<
        gbbs::uintE,
        gbbs::uintE,
        typename APSP<Graph>::Distance>>& expected) {
  using Distance = typename APSP<Graph>::Distance;

  CheckEqual(
      batch.size(),
      expected.size(),
      "distance-decrease batch size");

  for (size_t i = 0; i < expected.size(); ++i) {
    const auto& [u, v, d] = expected[i];

    CheckEqual(
        batch[i].p,
        u,
        "batch first endpoint");

    CheckEqual(
        batch[i].q,
        v,
        "batch second endpoint");

    CheckEqual(
        static_cast<Distance>(batch[i].new_distance),
        d,
        "batch new distance");
  }
}

template <class Graph>
void CommitBatch(
    APSP<Graph>& apsp,
    const gbbs::sequence<gbbs::kmedian::DistanceDecrease>& batch) {
  for (const auto& decrease : batch) {
    apsp.set_distance(
        decrease.p,
        decrease.q,
        decrease.new_distance);
  }
}

void TestUnweightedRandomSequence() {
  using W = gbbs::empty;
  using edge = std::tuple<gbbs::uintE, gbbs::uintE, W>;
  using Graph =
      gbbs::symmetric_graph<gbbs::symmetric_vertex, W>;

  constexpr gbbs::uintE n = 20;
  constexpr size_t kNumInsertions = 40;
  constexpr uint64_t kSeed = 42;

  std::mt19937_64 rng(kSeed);

  std::vector<edge> edges;
  std::vector<std::vector<bool>> present(
      n,
      std::vector<bool>(n, false));

  for (gbbs::uintE v = 1; v < n; ++v) {
    std::uniform_int_distribution<gbbs::uintE> parent_dist(0, v - 1);
    const gbbs::uintE u = parent_dist(rng);

    edges.emplace_back(u, v, W{});
    present[u][v] = true;
    present[v][u] = true;
  }

  auto MakeGraph = [&]() {
    gbbs::sequence<edge> seq(
        edges.begin(),
        edges.end());

    return Graph::from_edges(seq, n);
  };

  Graph G = MakeGraph();
  APSP<Graph> apsp(G);

  auto expected_current = FreshAPSP(G);
  VerifyMatrix(
      apsp,
      expected_current,
      "initial unweighted APSP");

  for (size_t update = 0; update < kNumInsertions; ++update) {
    gbbs::uintE u;
    gbbs::uintE v;

    do {
      u = static_cast<gbbs::uintE>(rng() % n);
      v = static_cast<gbbs::uintE>(rng() % n);
    } while (u == v || present[u][v]);

    if (u > v) {
      std::swap(u, v);
    }

    const auto before = expected_current;

    auto batch =
        apsp.GenerateBatch(
            G,
            u,
            v,
            W{});

    // GenerateBatch must leave the old metric untouched.
    VerifyMatrix(
        apsp,
        before,
        "transactional unweighted GenerateBatch");

    edges.emplace_back(u, v, W{});
    present[u][v] = true;
    present[v][u] = true;

    Graph updated = MakeGraph();
    auto after = FreshAPSP(updated);

    VerifyBatch<Graph>(
        batch,
        ExpectedDecreases(before, after));

    CommitBatch(
        apsp,
        batch);

    VerifyMatrix(
        apsp,
        after,
        "committed unweighted APSP");

    expected_current = std::move(after);
    G = std::move(updated);
  }

  std::cout
      << "IncrementalAPSP unweighted randomized insertion test passed."
      << std::endl;
}

void TestWeightedSequence() {
  using W = uint32_t;
  using edge = std::tuple<gbbs::uintE, gbbs::uintE, W>;
  using Graph =
      gbbs::symmetric_graph<gbbs::symmetric_vertex, W>;

  constexpr gbbs::uintE n = 8;

  std::vector<edge> edges = {
      edge{0, 1, 4},
      edge{1, 2, 3},
      edge{2, 3, 5},
      edge{3, 4, 2},
      edge{4, 5, 6},
      edge{5, 6, 3},
      edge{6, 7, 4},
      edge{0, 7, 30}
  };

  auto MakeGraph = [&]() {
    gbbs::sequence<edge> seq(
        edges.begin(),
        edges.end());

    return Graph::from_edges(seq, n);
  };

  Graph G = MakeGraph();
  APSP<Graph> apsp(G);

  const std::vector<std::tuple<gbbs::uintE, gbbs::uintE, W>>
      insertions = {
          {0, 4, 3},
          {2, 6, 2},
          {1, 7, 5},
          {3, 7, 1}
      };

  auto expected_current = FreshAPSP(G);

  for (const auto& [u, v, w] : insertions) {
    const auto before = expected_current;

    auto batch =
        apsp.GenerateBatch(
            G,
            u,
            v,
            w);

    VerifyMatrix(
        apsp,
        before,
        "transactional weighted GenerateBatch");

    edges.emplace_back(u, v, w);

    Graph updated = MakeGraph();
    auto after = FreshAPSP(updated);

    VerifyBatch<Graph>(
        batch,
        ExpectedDecreases(before, after));

    CommitBatch(
        apsp,
        batch);

    VerifyMatrix(
        apsp,
        after,
        "committed weighted APSP");

    expected_current = std::move(after);
    G = std::move(updated);
  }

  std::cout
      << "IncrementalAPSP weighted insertion test passed."
      << std::endl;
}

}  // namespace

int main() {
  TestUnweightedRandomSequence();
  TestWeightedSequence();

  return 0;
}

