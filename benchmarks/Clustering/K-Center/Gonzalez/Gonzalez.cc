// flags:
//   required:
//     -k <int>      : number of centers
//   optional:
//     -rounds <int> : number of times to run the algorithm
//     -s            : input graph is symmetric / undirected
//     -m            : input graph should be mmap'd
//     -c            : compressed graph input
//     -b            : binary graph input
//     -sc           : use single-core SSSP backend
//     -delta <dbl>  : bucket width 
//     -nb <int>     : number of buckets; must be a power of two
//     -seed <int>   : random seed for choosing the first center

#include <algorithm>
#include <ctime>
#include <type_traits>

#include "gbbs/gbbs.h"
#include "benchmarks/Clustering/K-Center/common/DispatchMain.h"
#include "Gonzalez.h"
#include "benchmarks/Clustering/K-Center/common/KCenterQuality.h"

namespace gbbs {

template <class Graph>
double Gonzalez_runner(Graph& G, commandLine P) {
  using W = typename Graph::weight_type;

  size_t k = P.getOptionLongValue("-k", 10);
  size_t num_buckets = P.getOptionLongValue("-nb", 32);
  double delta = P.getOptionDoubleValue("-delta", 1.0);
  uint64_t seed = static_cast<uint64_t>(
      P.getOptionLongValue("-seed", static_cast<long>(std::time(nullptr))));
  bool single_core = P.getOption("-sc");

  std::cout << "### Application: Gonzalez" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: "
            << parlay::num_workers()
            << std::endl;
  std::cout << "### Graph Type: "
            << (std::is_same<W, gbbs::empty>::value ? "unweighted" : "weighted")
            << std::endl;
  std::cout << "### SSSP Mode: "
            << (single_core ? "single-core" : "parallel") << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: -k = " << k
            << ", -delta = " << delta
            << ", -nb = " << num_buckets
            << ", -seed = " << seed
            << (single_core ? ", -sc = true" : ", -sc = false")
            << std::endl;
  std::cout << "### ------------------------------------" << std::endl;

  if constexpr (!std::is_same<W, gbbs::empty>::value) {
    if (!single_core &&
        num_buckets != (((uintE)1) << parlay::log2_up(num_buckets))) {
      std::cout << "Please specify a number of buckets that is a power of two\n";
      exit(-1);
    }
  }

  timer t;
  t.start();
  auto result = Gonzalez(G, k, delta, num_buckets, seed, single_core);
  double tt = t.stop();

  size_t preview = result.centers.size();
  if (preview > 0) {
    std::cout << "### Centers (first " << preview << "): ";
    for (size_t i = 0; i < preview; ++i) {
      std::cout << result.centers[i] << (i + 1 == preview ? "\n" : " ");
    }
  }

  auto q = gbbs::kcenter_quality::compute(
    G, result.centers, delta, num_buckets, single_core);

  std::cout << "num_centers = " << result.centers.size() << std::endl;
  gbbs::kcenter_quality::print(q);
  std::cout << "### Running Time: " << tt << std::endl;

  return tt;
}

}  // namespace gbbs

int main(int argc, char* argv[]) {
  auto app = [](auto& G, gbbs::commandLine P) {
    return gbbs::Gonzalez_runner(G, P);
  };
  return gbbs::kcenter_common::dispatch_main(argc, argv, app, false);
}