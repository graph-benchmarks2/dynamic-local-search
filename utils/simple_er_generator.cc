#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <string>
#include <vector>

namespace {

uint64_t splitmix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}


// -----------------------------------------------------------------------------
// Disjoint-set union for connectivity testing.
// -----------------------------------------------------------------------------

class DSU {
 public:
  explicit DSU(size_t n)
      : parent_(n),
        rank_(n, 0),
        components_(n) {
    std::iota(parent_.begin(), parent_.end(), size_t{0});
  }

  size_t Find(size_t x) {
    while (parent_[x] != x) {
      parent_[x] = parent_[parent_[x]];
      x = parent_[x];
    }
    return x;
  }

  void Union(size_t a, size_t b) {
    a = Find(a);
    b = Find(b);

    if (a == b) return;

    if (rank_[a] < rank_[b]) {
      std::swap(a, b);
    }

    parent_[b] = a;

    if (rank_[a] == rank_[b]) {
      ++rank_[a];
    }

    --components_;
  }

  bool Connected() const {
    return components_ == 1;
  }

 private:
  std::vector<size_t> parent_;
  std::vector<uint8_t> rank_;
  size_t components_;
};


void usage(const char* name) {
  std::cerr
      << "Usage:\n"
      << "  " << name
      << " <n> <m_undirected> <seed> <max_weight> "
      << "<out_unweighted> <out_weighted> [mode]\n\n"
      << "Modes:\n"
      << "  er            ordinary G(n,m) Erdős-Rényi graph (default)\n"
      << "  connected-er  G(n,m) conditioned on connectivity via rejection\n"
      << "                sampling\n"
      << "  tree          random spanning tree plus random additional edges\n\n"
      << "Writes exactly m_undirected distinct simple undirected edges u v,\n"
      << "and a weighted copy u v w with weights in [1,max_weight].\n";
}


// -----------------------------------------------------------------------------
// Ordinary G(n,m) edge generation.
//
// This is essentially the original generator.
// -----------------------------------------------------------------------------

std::vector<uint64_t> GenerateRandomEdges(
    uint64_t n,
    uint64_t m,
    std::mt19937_64& rng,
    bool print_progress) {

  std::uniform_int_distribution<uint64_t> vertex_dist(0, n - 1);

  std::vector<uint64_t> edges;
  edges.reserve(m);

  while (edges.size() < m) {
    const uint64_t need = m - edges.size();

    const uint64_t batch =
        std::max<uint64_t>(
            need + need / 100 + 1024,
            4096);

    for (uint64_t i = 0; i < batch; ++i) {
      uint64_t u = vertex_dist(rng);
      uint64_t v = vertex_dist(rng);

      if (u == v) continue;

      if (u > v) {
        std::swap(u, v);
      }

      edges.push_back(u * n + v);
    }

    std::sort(edges.begin(), edges.end());

    edges.erase(
        std::unique(edges.begin(), edges.end()),
        edges.end());

    if (print_progress) {
      std::cerr << "distinct_edges="
                << edges.size()
                << " target="
                << m
                << std::endl;
    }
  }

  if (edges.size() > m) {
    std::shuffle(edges.begin(), edges.end(), rng);
    edges.resize(m);
    std::sort(edges.begin(), edges.end());
  }

  return edges;
}


// -----------------------------------------------------------------------------
// Connectivity check directly on encoded edge list.
// -----------------------------------------------------------------------------

bool IsConnected(
    uint64_t n,
    const std::vector<uint64_t>& edges) {

  DSU dsu(n);

  for (uint64_t key : edges) {
    const uint64_t u = key / n;
    const uint64_t v = key % n;

    dsu.Union(
        static_cast<size_t>(u),
        static_cast<size_t>(v));
  }

  return dsu.Connected();
}


// -----------------------------------------------------------------------------
// Uniform random labeled spanning tree using a random Prüfer sequence.
//
// Every labeled tree on {0,...,n-1} occurs with equal probability.
// -----------------------------------------------------------------------------

std::vector<uint64_t> GenerateUniformRandomTree(
    uint64_t n,
    std::mt19937_64& rng) {

  std::vector<uint64_t> edges;
  edges.reserve(n - 1);

  if (n == 2) {
    edges.push_back(1);
    return edges;
  }

  std::uniform_int_distribution<uint64_t> vertex_dist(0, n - 1);

  // Prüfer sequence has length n-2.
  std::vector<uint64_t> prufer(n - 2);

  std::vector<uint64_t> degree(n, 1);

  for (uint64_t i = 0; i < n - 2; ++i) {
    prufer[i] = vertex_dist(rng);
    ++degree[prufer[i]];
  }

  std::priority_queue<
      uint64_t,
      std::vector<uint64_t>,
      std::greater<uint64_t>>
      leaves;

  for (uint64_t v = 0; v < n; ++v) {
    if (degree[v] == 1) {
      leaves.push(v);
    }
  }

  for (uint64_t p : prufer) {
    const uint64_t leaf = leaves.top();
    leaves.pop();

    uint64_t u = leaf;
    uint64_t v = p;

    if (u > v) {
      std::swap(u, v);
    }

    edges.push_back(u * n + v);

    --degree[leaf];
    --degree[p];

    if (degree[p] == 1) {
      leaves.push(p);
    }
  }

  uint64_t u = leaves.top();
  leaves.pop();

  uint64_t v = leaves.top();
  leaves.pop();

  if (u > v) {
    std::swap(u, v);
  }

  edges.push_back(u * n + v);

  std::sort(edges.begin(), edges.end());

  return edges;
}


// -----------------------------------------------------------------------------
// Connected sparse graph:
//
//   1. generate a uniform random labeled spanning tree,
//   2. add random distinct edges until exactly m edges exist.
//
// This is NOT G(n,m) conditioned on connectivity.
// -----------------------------------------------------------------------------

}  // namespace


int main(int argc, char** argv) {
  if (argc != 7 && argc != 8) {
    usage(argv[0]);
    return 1;
  }

  const uint64_t n =
      std::stoull(argv[1]);

  const uint64_t m =
      std::stoull(argv[2]);

  const uint64_t seed =
      std::stoull(argv[3]);

  const uint64_t max_weight =
      std::stoull(argv[4]);

  const std::string out_unweighted =
      argv[5];

  const std::string out_weighted =
      argv[6];

  const std::string mode =
      (argc == 8)
          ? argv[7]
          : "er";

  if (n < 2 || m == 0 || max_weight == 0) {
    std::cerr << "Invalid parameters.\n";
    return 1;
  }

  const long double max_edges =
      (static_cast<long double>(n) *
       static_cast<long double>(n - 1)) /
      2.0L;

  if (static_cast<long double>(m) > max_edges) {
    std::cerr
        << "Requested more edges than possible "
           "simple undirected edges.\n";
    return 1;
  }

  if (mode != "er" &&
      mode != "connected-er" &&
      mode != "tree") {

    std::cerr
        << "Unknown mode: "
        << mode
        << "\n";

    usage(argv[0]);
    return 1;
  }

  if ((mode == "connected-er" || mode == "tree") &&
      m < n - 1) {

    std::cerr
        << "A connected simple graph on n vertices "
           "requires at least n-1 edges.\n";
    return 1;
  }

  std::mt19937_64 rng(seed);

  std::vector<uint64_t> edges;


  // ---------------------------------------------------------------------------
  // Ordinary ER
  // ---------------------------------------------------------------------------

  if (mode == "er") {
    edges =
        GenerateRandomEdges(
            n,
            m,
            rng,
            true);
  }


  // ---------------------------------------------------------------------------
  // ER conditioned on connectivity
  // ---------------------------------------------------------------------------

  else if (mode == "connected-er") {
    uint64_t attempt = 0;

    while (true) {
      ++attempt;

      std::cerr
          << "connected-er attempt="
          << attempt
          << std::endl;

      edges =
          GenerateRandomEdges(
              n,
              m,
              rng,
              false);

      if (IsConnected(n, edges)) {
        std::cerr
            << "connected after "
            << attempt
            << " attempt(s)"
            << std::endl;

        break;
      }

      std::cerr
          << "rejected: graph disconnected"
          << std::endl;
    }
  }


  // ---------------------------------------------------------------------------
  // Tree + random additional edges
  // ---------------------------------------------------------------------------

  else {
    // Generate the tree separately so we can guarantee that these
    // n-1 edges are never removed.
    auto tree =
        GenerateUniformRandomTree(n, rng);

    edges = tree;

    std::uniform_int_distribution<uint64_t>
        vertex_dist(0, n - 1);

    while (edges.size() < m) {
      const uint64_t need =
          m - edges.size();

      const uint64_t batch =
          std::max<uint64_t>(
              need + need / 100 + 1024,
              4096);

      for (uint64_t i = 0; i < batch; ++i) {
        uint64_t u =
            vertex_dist(rng);

        uint64_t v =
            vertex_dist(rng);

        if (u == v) continue;

        if (u > v) {
          std::swap(u, v);
        }

        edges.push_back(u * n + v);
      }

      std::sort(edges.begin(), edges.end());

      edges.erase(
          std::unique(edges.begin(), edges.end()),
          edges.end());

      std::cerr
          << "distinct_edges="
          << edges.size()
          << " target="
          << m
          << std::endl;
    }

    /*
     * Because the batching above may overshoot m, we cannot randomly
     * truncate the entire edge set: that might delete a tree edge.
     *
     * Reconstruct a tree-membership marker, keep all tree edges, and
     * sample only enough of the remaining edges.
     */
    if (edges.size() > m) {
      std::vector<uint64_t> extra_edges;

      extra_edges.reserve(
          edges.size() - tree.size());

      for (uint64_t key : edges) {
        if (!std::binary_search(
                tree.begin(),
                tree.end(),
                key)) {

          extra_edges.push_back(key);
        }
      }

      std::shuffle(
          extra_edges.begin(),
          extra_edges.end(),
          rng);

      const uint64_t needed_extra =
          m - (n - 1);

      extra_edges.resize(needed_extra);

      edges = tree;

      edges.insert(
          edges.end(),
          extra_edges.begin(),
          extra_edges.end());

      std::sort(
          edges.begin(),
          edges.end());
    }
  }


  // ---------------------------------------------------------------------------
  // Sanity checks
  // ---------------------------------------------------------------------------

  if (edges.size() != m) {
    std::cerr
        << "Internal error: generated "
        << edges.size()
        << " edges instead of "
        << m
        << ".\n";

    return 1;
  }

  if (mode != "er" &&
      !IsConnected(n, edges)) {

    std::cerr
        << "Internal error: connected mode "
           "produced disconnected graph.\n";

    return 1;
  }


  // ---------------------------------------------------------------------------
  // Output
  // ---------------------------------------------------------------------------

  std::ofstream fout_u(out_unweighted);
  std::ofstream fout_w(out_weighted);

  if (!fout_u) {
    std::cerr
        << "Could not open unweighted output: "
        << out_unweighted
        << "\n";

    return 1;
  }

  if (!fout_w) {
    std::cerr
        << "Could not open weighted output: "
        << out_weighted
        << "\n";

    return 1;
  }

  for (uint64_t key : edges) {
    const uint64_t u =
        key / n;

    const uint64_t v =
        key % n;

    const uint64_t w =
        (splitmix64(
             key ^
             seed ^
             0xdeadbeefcafebabeULL) %
         max_weight) +
        1;

    fout_u
        << u << ' '
        << v << '\n';

    fout_w
        << u << ' '
        << v << ' '
        << w << '\n';
  }

  std::cerr
      << "done mode="
      << mode
      << " edges="
      << m
      << "\n";

  return 0;
}
