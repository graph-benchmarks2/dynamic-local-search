// Static k-median local-search launcher.
//
// Timing convention:
//   initializer SSSP time is reported separately because some initializers
//   still run graph SSSPs. Local search itself is metric-backed and therefore
//   reports only wall time. APSP preprocessing is metric infrastructure and is
//   excluded from the paper-facing algorithm time.

// flags:
//   required:
//     -k <int>        : number of centers
//   optional:
//     -rounds <int>   : number of times to run the algorithm
//     -s              : input graph is symmetric / undirected
//     -m              : input graph should be mmap'd
//     -c              : compressed graph input
//     -b              : binary graph input
//     -sc             : use single-core SSSP backend
//     -delta <dbl>    : DeltaStepping bucket width
//     -y <int>        : incoming-center candidate for swap-column debug test
//     -ls-delta <dbl> : local-search significance parameter (default: 0.1)
//     -nb <int>       : DeltaStepping number of buckets; power of two
//     -seed <int>     : seed for random center initialization
//     -init <string>  : random | greedy | gonzalez (default: random)
//     -verbose        : print per-vertex solution state and LS iterations
//     -sssp-verbose   : print diagnostics from every SSSP call
//     -debug-swaps    : run expensive swap-table correctness checks

#include <algorithm>
#include <ctime>
#include <type_traits>

#include "gbbs/gbbs.h"
#include "benchmarks/Clustering/K-Center/common/DispatchMain.h"
#include "StaticLocalSearch.h"
#include "benchmarks/Clustering/K-Median/common/IncrementalAPSP.h"
#include "benchmarks/Clustering/K-Median/common/Initializers.h"

namespace gbbs {

    template <class Graph>
    double StaticLocalSearch_runner(Graph& G, commandLine P) {
        using W = typename Graph::weight_type;

        size_t k = P.getOptionLongValue("-k", 10);
        size_t num_buckets = P.getOptionLongValue("-nb", 32);
        double delta = P.getOptionDoubleValue("-delta", 1.0);
        double ls_delta = P.getOptionDoubleValue("-ls-delta", 0.1);
        uint64_t seed = static_cast<uint64_t>(P.getOptionLongValue("-seed", 42));
        std::string initializer = P.getOptionValue("-init", "random");
        bool single_core = P.getOption("-sc");
        bool verbose = P.getOption("-verbose");
        bool sssp_verbose = P.getOption("-sssp-verbose");
        bool debug_swaps = P.getOption("-debug-swaps");

        std::cout << "### Application: StaticLocalSearch" << std::endl;
        std::cout << "### Graph: " << P.getArgument(0) << std::endl;
        std::cout << "### Threads: " << parlay::num_workers() << std::endl;
        std::cout << "### Graph Type: "
                  << (std::is_same<W, gbbs::empty>::value ? "unweighted" : "weighted")
                  << std::endl;
        std::cout << "### SSSP Mode: " << (single_core ? "single-core" : "parallel") << std::endl;
        std::cout << "### n: " << G.n << std::endl;
        std::cout << "### m: " << G.m << std::endl;
        std::cout << "### Params: -k = " << k
                  << ", -delta = " << delta
                  << ", -ls-delta = " << ls_delta
                  << ", -nb = " << num_buckets
                  << ", -seed = " << seed
                  << ", -init = " << initializer
                  << (single_core ? ", -sc = true" : ", -sc = false")
                  << std::endl;
        std::cout << "### ------------------------------------" << std::endl;

        if (k == 0) {
            std::cout << "ERROR: k must be at least 1." << std::endl;
            std::exit(-1);
        }

        if (k > G.n) {
            std::cout << "ERROR: k = " << k << " exceeds n = " << G.n << "." << std::endl;
            std::exit(-1);
        }

        if constexpr (!std::is_same<W, gbbs::empty>::value) {
            if (!single_core && num_buckets != (((uintE)1) << parlay::log2_up(num_buckets))) {
                std::cout << "Please specify a number of buckets that is a power of two\n";
                std::exit(-1);
            }
        }

        kmedian::SSSPTiming initializer_sssp_timing;
        timer initializer_timer;
        initializer_timer.start();

        auto centers = kmedian::InitializeCenters(
            G, initializer, k, delta, num_buckets, seed, single_core, sssp_verbose, &initializer_sssp_timing);

        const double initializer_wall_time = initializer_timer.stop();
        const double initializer_distance_time = initializer_sssp_timing.total_time;
        const double initializer_clustering_time =
            std::max(0.0, initializer_wall_time - initializer_distance_time);

        // Exact APSP preprocessing.
        kmedian::SSSPTiming apsp_sssp_timing;
        timer apsp_timer;
        apsp_timer.start();

        kmedian::IncrementalAPSP<Graph> apsp(
            G,
            std::numeric_limits<double>::infinity(),
            &apsp_sssp_timing);

        const double apsp_wall_time =
            apsp_timer.stop();

        // State construction from the already available metric.
        timer state_timer;
        state_timer.start();

        auto state =
            kmedian::BuildSolutionStateFromDistances(
                apsp,
                centers);

        const double state_wall_time =
            state_timer.stop();

        std::cout << "### Centers: ";
        for (size_t i = 0; i < state.centers.size(); ++i) {
            std::cout << state.centers[i] << (i + 1 == state.centers.size() ? "\n" : " ");
        }

        std::cout << "initializer = " << initializer << std::endl;
        std::cout << "initializer_wall_time = " << initializer_wall_time << std::endl;
        std::cout << "initializer_distance_time = " << initializer_distance_time << std::endl;
        std::cout << "initializer_clustering_time = " << initializer_clustering_time << std::endl;
        std::cout << "apsp_preprocessing_wall_time = " << apsp_wall_time << std::endl;
        std::cout << "apsp_preprocessing_sssp_time = " << apsp_sssp_timing.total_time << std::endl;
        std::cout << "apsp_preprocessing_sssp_calls = " << apsp_sssp_timing.calls << std::endl;
        std::cout << "solution_state_wall_time = " << state_wall_time << std::endl;
        std::cout << "cost = " << state.cost << std::endl;
        std::cout << "unreachable_vertices = " << state.unreachable_vertices << std::endl;

        if (verbose) {
            std::cout << "### Vertex state:" << std::endl;

            for (size_t v = 0; v < G.n; ++v) {
                std::cout << "v=" << v
                          << " nearest_center=" << state.nearest_center[v]
                          << " nearest_dist=" << state.nearest_dist[v]
                          << " second_nearest_center=" << state.second_nearest_center[v]
                          << " second_nearest_dist=" << state.second_nearest_dist[v]
                          << std::endl;
            }
        }

        if (debug_swaps) {
            uintE y = kmedian::kNoCenter;
            long requested_y = P.getOptionLongValue("-y", -1);

            if (requested_y >= 0) {
                y = static_cast<uintE>(requested_y);

                if (y >= G.n) {
                    std::cout << "ERROR: -y " << y << " is outside the graph." << std::endl;
                    std::exit(-1);
                }

                if (state.is_center[y]) {
                    std::cout << "ERROR: -y " << y << " is already a center." << std::endl;
                    std::exit(-1);
                }
            } else {
                // Default: smallest vertex ID that is not currently a center.
                for (uintE v = 0; v < G.n; ++v) {
                    if (!state.is_center[v]) {
                        y = v;
                        break;
                    }
                }
            }

            if (y == kmedian::kNoCenter) {
                std::cout << "ERROR: No non-center exists." << std::endl;
                std::exit(-1);
            }

            auto column =
                kmedian::ComputeSwapColumnFromDistances(
                    apsp,
                    state,
                    y);

            std::cout << "### Swap table column y=" << y << std::endl;

            for (size_t row = 0; row < state.centers.size(); ++row) {
                std::cout << "remove=" << state.centers[row]
                          << " add=" << y
                          << " gain=" << column.gains[row]
                          << std::endl;
            }

            std::cout << "best_outgoing_center = " << column.best_outgoing_center << std::endl;
            std::cout << "column_max_gain = " << column.max_gain << std::endl;

            kmedian::VerifySwapColumnBruteForce(
                apsp,
                state,
                column);

            std::cout << "### ====================================" << std::endl;
            std::cout << "### Full static swap table" << std::endl;

            kmedian::StaticSwapTable table(
                apsp,
                state);

            std::cout << "### Rows (outgoing centers): ";
            for (size_t row = 0; row < table.num_rows(); ++row) {
                std::cout << table.outgoing_center(row)
                          << (row + 1 == table.num_rows() ? "\n" : " ");
            }

            std::cout << "### Columns (incoming centers): ";
            for (size_t col = 0; col < table.num_cols(); ++col) {
                std::cout << table.incoming_center(col)
                          << (col + 1 == table.num_cols() ? "\n" : " ");
            }

            for (size_t row = 0; row < table.num_rows(); ++row) {
                const uintE x = table.outgoing_center(row);

                for (size_t col = 0; col < table.num_cols(); ++col) {
                    const uintE incoming = table.incoming_center(col);
                    std::cout << "Delta(" << x << "," << incoming << ") = "
                              << table.gain(row, col)
                              << std::endl;
                }
            }

            std::cout << "global_best_remove = " << table.best_outgoing_center() << std::endl;
            std::cout << "global_best_add = " << table.best_incoming_center() << std::endl;
            std::cout << "global_max_gain = " << table.max_gain() << std::endl;

            kmedian::VerifyStaticSwapTableBruteForce(
                apsp,
                state,
                table);
        }

        std::cout << "### ====================================" << std::endl;
        std::cout << "### Static local search" << std::endl;

        auto ls_result =
            kmedian::RunStaticLocalSearchFromDistances(
                apsp,
                centers,
                ls_delta,
                verbose);

        std::cout << "### Static local-search result" << std::endl;
        std::cout << "initial_cost = " << ls_result.initial_cost << std::endl;
        std::cout << "final_cost = " << ls_result.final_cost << std::endl;
        std::cout << "num_swaps = " << ls_result.num_swaps << std::endl;
        std::cout << "final_max_gain = " << ls_result.final_max_gain << std::endl;
        std::cout << "final_threshold = " << ls_result.final_threshold << std::endl;
        std::cout << "final_centers = ";

        for (size_t i = 0; i < ls_result.state.centers.size(); ++i) {
            std::cout << ls_result.state.centers[i]
                      << (i + 1 == ls_result.state.centers.size() ? "\n" : " ");
        }

        std::cout << "local_search_wall_time = " << ls_result.wall_time << std::endl;

        const double total_algorithm_wall_time =
            initializer_wall_time + ls_result.wall_time;
        const double total_algorithm_distance_time =
            initializer_distance_time;
        const double total_algorithm_clustering_time =
            initializer_clustering_time + ls_result.wall_time;

        std::cout << "total_algorithm_wall_time = " << total_algorithm_wall_time << std::endl;
        std::cout << "total_algorithm_distance_time = " << total_algorithm_distance_time << std::endl;
        std::cout << "total_algorithm_clustering_time = " << total_algorithm_clustering_time << std::endl;

        return total_algorithm_clustering_time;
    }

}  // namespace gbbs

int main(int argc, char* argv[]) {
    auto app = [](auto& G, gbbs::commandLine P) {
        return gbbs::StaticLocalSearch_runner(G, P);
    };

    return gbbs::kcenter_common::dispatch_main(argc, argv, app, false);
}
