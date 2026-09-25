#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "gbbs/gbbs.h"

#include "benchmarks/Clustering/K-Center/common/DispatchMain.h"
#include "benchmarks/Clustering/K-Median/common/IncrementalAPSP.h"
#include "benchmarks/Clustering/K-Median/common/Initializers.h"
#include "benchmarks/Clustering/K-Median/StaticLocalSearch/StaticLocalSearch.h"
#include "benchmarks/Clustering/K-Median/StaticLocalSearch/WarmStartDriver.h"

namespace gbbs {

    template <class Graph>
    double WarmStartStaticLocalSearch_runner(Graph& G, commandLine P) {
        const size_t k = P.getOptionLongValue("-k", 10);
        const double delta = P.getOptionDoubleValue("-delta", 1.0);
        const double ls_delta = P.getOptionDoubleValue("-ls-delta", 0.1);
        const size_t num_buckets = P.getOptionLongValue("-nb", 32);
        const uint64_t seed = static_cast<uint64_t>(P.getOptionLongValue("-seed", 42));
        const std::string initializer = P.getOptionValue("-init", "random");
        const std::string fixed_centers_string = P.getOptionValue("-centers", "");
        const std::string updates_file = P.getOptionValue("-updates", "");
        const bool single_core = P.getOption("-sc");
        const bool verbose = P.getOption("-verbose");
        const bool sssp_verbose = P.getOption("-sssp-verbose");
        const size_t cold_quality_trials =
            P.getOptionLongValue("-quality-cold-trials", 0);
        const uint64_t cold_quality_seed = static_cast<uint64_t>(
            P.getOptionLongValue("-quality-cold-seed", seed));
        const std::string cold_quality_checkpoints_string =
            P.getOptionValue("-quality-cold-checkpoints", "");
        std::vector<size_t> cold_quality_checkpoints;

        if (!cold_quality_checkpoints_string.empty()) {
            std::stringstream ss(cold_quality_checkpoints_string);
            std::string token;

            while (std::getline(ss, token, ',')) {
                if (token.empty()) {
                    std::cout << "ERROR: Invalid -quality-cold-checkpoints specification: "
                              << cold_quality_checkpoints_string << std::endl;
                    std::exit(-1);
                }

                uint64_t value;
                try {
                    value = std::stoull(token);
                } catch (...) {
                    std::cout << "ERROR: Invalid quality checkpoint: " << token << std::endl;
                    std::exit(-1);
                }

                if (value == 0) {
                    std::cout << "ERROR: Quality checkpoints are 1-based and must be positive."
                              << std::endl;
                    std::exit(-1);
                }

                cold_quality_checkpoints.push_back(static_cast<size_t>(value));
            }

            std::sort(cold_quality_checkpoints.begin(), cold_quality_checkpoints.end());
            cold_quality_checkpoints.erase(
                std::unique(
                    cold_quality_checkpoints.begin(),
                    cold_quality_checkpoints.end()),
                cold_quality_checkpoints.end());
        }

        if (updates_file.empty()) {
            std::cout << "ERROR: Warm-start driver requires -updates <file>." << std::endl;
            std::exit(-1);
        }

        if (k == 0) {
            std::cout << "ERROR: k must be positive." << std::endl;
            std::exit(-1);
        }

        if (k > G.n) {
            std::cout << "ERROR: k = " << k << " exceeds n = " << G.n << "." << std::endl;
            std::exit(-1);
        }

        parlay::sequence<uintE> fixed_centers;

        if (!fixed_centers_string.empty()) {
            std::stringstream ss(fixed_centers_string);
            std::string token;
            std::unordered_set<uintE> seen;

            while (std::getline(ss, token, ',')) {
                if (token.empty()) {
                    std::cout << "ERROR: Invalid -centers specification: "
                              << fixed_centers_string
                              << std::endl;
                    std::exit(-1);
                }

                uint64_t value;

                try {
                    value = std::stoull(token);
                } catch (...) {
                    std::cout << "ERROR: Invalid center vertex: " << token << std::endl;
                    std::exit(-1);
                }

                if (value >= G.n) {
                    std::cout << "ERROR: Center " << value
                              << " is outside vertex range [0, " << G.n - 1 << "]."
                              << std::endl;
                    std::exit(-1);
                }

                const uintE center = static_cast<uintE>(value);
                if (!seen.insert(center).second) {
                    std::cout << "ERROR: Duplicate center " << center << " in -centers." << std::endl;
                    std::exit(-1);
                }

                fixed_centers.push_back(center);
            }

            if (fixed_centers.size() != k) {
                std::cout << "ERROR: -centers contains " << fixed_centers.size()
                          << " centers, but k = " << k << "."
                          << std::endl;
                std::exit(-1);
            }
        }

        std::cout << "### Application: WarmStartStaticLocalSearch" << std::endl;
        std::cout << "### n: " << G.n << std::endl;
        std::cout << "### m: " << G.m << std::endl;

        if (!fixed_centers.empty()) {
            std::cout << "### initializer = fixed" << std::endl;
            std::cout << "### fixed centers =";

            for (const auto center : fixed_centers) {
                std::cout << " " << center;
            }

            std::cout << std::endl;
        } else {
            std::cout << "### initializer = " << initializer << std::endl;
        }

        std::cout << "### updates = " << updates_file << std::endl;
        std::cout << "### ls_delta = " << ls_delta << std::endl;

        if (cold_quality_trials > 0) {
            std::cout << "### quality_cold_trials = " << cold_quality_trials << std::endl;
            std::cout << "### quality_cold_seed = " << cold_quality_seed << std::endl;
            std::cout << "### quality baseline = fresh random cold starts" << std::endl;
            if (!cold_quality_checkpoints.empty()) {
                std::cout << "### quality_cold_checkpoints =";
                for (const auto checkpoint : cold_quality_checkpoints) {
                    std::cout << " " << checkpoint;
                }
                std::cout << std::endl;
            }
        }

        std::cout << "### ------------------------------------" << std::endl;

        parlay::sequence<uintE> initial_centers;
        kmedian::SSSPTiming initializer_sssp_timing;

        timer initializer_timer;
        initializer_timer.start();

        if (!fixed_centers.empty()) {
            initial_centers = fixed_centers;
        } else {
            initial_centers = kmedian::InitializeCenters(
                G,
                initializer,
                k,
                delta,
                num_buckets,
                seed,
                single_core,
                sssp_verbose,
                &initializer_sssp_timing);
        }

        const double initializer_wall_time = initializer_timer.stop();
        const double initializer_distance_time = initializer_sssp_timing.total_time;
        const double initializer_clustering_time =
            std::max(0.0, initializer_wall_time - initializer_distance_time);

        using W = typename Graph::weight_type;
        using edge = std::tuple<uintE, uintE, W>;

        auto base_edges = kmedian::ExtractUniqueUndirectedEdges(G);
        parlay::sequence<edge> apsp_edges(base_edges.begin(), base_edges.end());
        auto apsp_graph = kmedian::BuildGraphFromEdges<W>(apsp_edges, G.n);

        kmedian::SSSPTiming apsp_sssp_timing;

        timer apsp_timer;
        apsp_timer.start();

        kmedian::IncrementalAPSP<decltype(apsp_graph)> incremental_apsp(
            apsp_graph,
            std::numeric_limits<double>::infinity(),
            &apsp_sssp_timing);

        const double initial_apsp_wall_time = apsp_timer.stop();

        auto initial_ls = kmedian::RunStaticLocalSearchFromDistances(
            incremental_apsp,
            initial_centers,
            ls_delta,
            verbose);

        const double initial_solution_wall_time =
            initializer_wall_time + initial_apsp_wall_time + initial_ls.wall_time;
        const double initial_solution_distance_time =
            initializer_distance_time + apsp_sssp_timing.total_time;
        const double initial_solution_clustering_time =
            initializer_clustering_time + initial_ls.wall_time;

        std::cout << "### Initial solution" << std::endl;
        std::cout << "initial_final_cost = " << initial_ls.final_cost << std::endl;
        std::cout << "initial_num_swaps = " << initial_ls.num_swaps << std::endl;
        std::cout << "initial_initializer_wall_time = " << initializer_wall_time << std::endl;
        std::cout << "initial_initializer_distance_time = " << initializer_distance_time << std::endl;
        std::cout << "initial_initializer_clustering_time = " << initializer_clustering_time << std::endl;
        std::cout << "initial_apsp_wall_time = " << initial_apsp_wall_time << std::endl;
        std::cout << "initial_apsp_sssp_time = " << apsp_sssp_timing.total_time << std::endl;
        std::cout << "initial_apsp_sssp_calls = " << apsp_sssp_timing.calls << std::endl;
        std::cout << "initial_local_search_wall_time = " << initial_ls.wall_time << std::endl;
        std::cout << "initial_solution_wall_time = " << initial_solution_wall_time << std::endl;
        std::cout << "initial_solution_distance_time = " << initial_solution_distance_time << std::endl;
        std::cout << "initial_solution_clustering_time = " << initial_solution_clustering_time << std::endl;
        std::cout << "### ====================================" << std::endl;

        timer update_timer;
        update_timer.start();

        auto update_results = kmedian::RunWarmStartSequence(
            apsp_graph,
            incremental_apsp,
            initial_ls.state.centers,
            initial_ls.final_cost,
            updates_file,
            ls_delta,
            verbose,
            cold_quality_trials,
            cold_quality_seed,
            cold_quality_checkpoints);

        const double driver_update_wall_time = update_timer.stop();

        size_t total_swaps = 0;
        double total_update_wall_time = 0.0;
        size_t total_swap_table_builds = 0;
        double total_initial_state_build_wall_time = 0.0;
        double total_swap_table_build_wall_time = 0.0;
        double total_state_rebuild_after_swap_wall_time = 0.0;
        double total_cold_initializer_wall_time = 0.0;
        double total_cold_local_search_wall_time = 0.0;
        double total_quality_baseline_wall_time = 0.0;

        for (const auto& r : update_results) {
            total_swaps += r.num_swaps;
            total_update_wall_time += r.update_wall_time;
            total_swap_table_builds += r.num_swap_table_builds;
            total_initial_state_build_wall_time +=
                r.initial_state_build_wall_time;
            total_swap_table_build_wall_time +=
                r.swap_table_build_wall_time;
            total_state_rebuild_after_swap_wall_time +=
                r.state_rebuild_after_swap_wall_time;
            total_cold_initializer_wall_time += r.cold_initializer_wall_time;
            total_cold_local_search_wall_time += r.cold_local_search_wall_time;
            total_quality_baseline_wall_time += r.quality_baseline_wall_time;
        }

        std::cout << "### Warm-start summary" << std::endl;
        std::cout << "num_updates = " << update_results.size() << std::endl;
        std::cout << "total_update_swaps = " << total_swaps << std::endl;
        std::cout << "total_update_wall_time = " << total_update_wall_time << std::endl;
        std::cout << "total_swap_table_builds = "
                  << total_swap_table_builds << std::endl;
        std::cout << "total_initial_state_build_wall_time = "
                  << total_initial_state_build_wall_time << std::endl;
        std::cout << "total_swap_table_build_wall_time = "
                  << total_swap_table_build_wall_time << std::endl;
        std::cout << "total_state_rebuild_after_swap_wall_time = "
                  << total_state_rebuild_after_swap_wall_time << std::endl;
        if (cold_quality_trials == 0) {
            std::cout << "driver_update_wall_time_including_apsp_and_graph_rebuilds = "
                      << driver_update_wall_time
                      << std::endl;
        } else {
            std::cout << "driver_wall_time_including_quality_baselines = "
                      << driver_update_wall_time
                      << std::endl;
            std::cout << "driver_operational_wall_time_excluding_quality_baselines = "
                      << std::max(
                             0.0,
                             driver_update_wall_time - total_quality_baseline_wall_time)
                      << std::endl;
        }

        if (cold_quality_trials > 0) {
            double same_seed_ratio_sum = 0.0;
            double best_ratio_sum = 0.0;
            double max_same_seed_ratio = 0.0;
            double max_best_ratio = 0.0;
            size_t maintained_better_same_seed = 0;
            size_t maintained_equal_same_seed = 0;
            size_t maintained_worse_same_seed = 0;
            size_t maintained_better_best = 0;
            size_t maintained_equal_best = 0;
            size_t maintained_worse_best = 0;
            std::vector<double> same_seed_ratios;
            std::vector<double> best_ratios;
            size_t quality_evaluations = 0;
            same_seed_ratios.reserve(update_results.size());
            best_ratios.reserve(update_results.size());

            const auto classify = [](double maintained,
                                     double cold,
                                     size_t& better,
                                     size_t& equal,
                                     size_t& worse) {
                constexpr double kTolerance = 1e-9;
                if (maintained < cold - kTolerance) {
                    ++better;
                } else if (maintained > cold + kTolerance) {
                    ++worse;
                } else {
                    ++equal;
                }
            };

            for (const auto& r : update_results) {
                if (!r.quality_evaluated) {
                    continue;
                }

                ++quality_evaluations;
                same_seed_ratio_sum += r.maintained_over_cold_same_seed;
                best_ratio_sum += r.maintained_over_cold_best;
                max_same_seed_ratio =
                    std::max(max_same_seed_ratio, r.maintained_over_cold_same_seed);
                max_best_ratio =
                    std::max(max_best_ratio, r.maintained_over_cold_best);
                same_seed_ratios.push_back(r.maintained_over_cold_same_seed);
                best_ratios.push_back(r.maintained_over_cold_best);

                classify(
                    r.final_cost,
                    r.cold_same_seed_cost,
                    maintained_better_same_seed,
                    maintained_equal_same_seed,
                    maintained_worse_same_seed);
                classify(
                    r.final_cost,
                    r.cold_best_cost,
                    maintained_better_best,
                    maintained_equal_best,
                    maintained_worse_best);
            }

            std::sort(same_seed_ratios.begin(), same_seed_ratios.end());
            std::sort(best_ratios.begin(), best_ratios.end());

            const auto median = [](const std::vector<double>& values) {
                if (values.empty()) {
                    return 0.0;
                }

                const size_t mid = values.size() / 2;
                if (values.size() % 2 == 1) {
                    return values[mid];
                }

                return 0.5 * (values[mid - 1] + values[mid]);
            };

            const double count = static_cast<double>(quality_evaluations);

            std::cout << "### Quality summary" << std::endl;
            std::cout << "quality_cold_trials = " << cold_quality_trials << std::endl;
            std::cout << "quality_cold_seed = " << cold_quality_seed << std::endl;
            std::cout << "quality_evaluated_checkpoints = " << quality_evaluations << std::endl;
            std::cout << "mean_maintained_over_cold_same_seed = "
                      << same_seed_ratio_sum / count << std::endl;
            std::cout << "median_maintained_over_cold_same_seed = "
                      << median(same_seed_ratios) << std::endl;
            std::cout << "max_maintained_over_cold_same_seed = "
                      << max_same_seed_ratio << std::endl;
            std::cout << "maintained_better_than_cold_same_seed = "
                      << maintained_better_same_seed << std::endl;
            std::cout << "maintained_equal_to_cold_same_seed = "
                      << maintained_equal_same_seed << std::endl;
            std::cout << "maintained_worse_than_cold_same_seed = "
                      << maintained_worse_same_seed << std::endl;
            std::cout << "mean_maintained_over_cold_best = "
                      << best_ratio_sum / count << std::endl;
            std::cout << "median_maintained_over_cold_best = "
                      << median(best_ratios) << std::endl;
            std::cout << "max_maintained_over_cold_best = "
                      << max_best_ratio << std::endl;
            std::cout << "maintained_better_than_cold_best = "
                      << maintained_better_best << std::endl;
            std::cout << "maintained_equal_to_cold_best = "
                      << maintained_equal_best << std::endl;
            std::cout << "maintained_worse_than_cold_best = "
                      << maintained_worse_best << std::endl;
            std::cout << "cold_initializer_wall_time = "
                      << total_cold_initializer_wall_time << std::endl;
            std::cout << "cold_local_search_wall_time = "
                      << total_cold_local_search_wall_time << std::endl;
            std::cout << "cold_algorithm_wall_time = "
                      << total_cold_initializer_wall_time +
                             total_cold_local_search_wall_time
                      << std::endl;
            std::cout << "quality_baseline_wall_time = "
                      << total_quality_baseline_wall_time << std::endl;
        }

        // Warm local search's wall time is also the clustering update time.
        return total_update_wall_time;
    }

}  // namespace gbbs

int main(int argc, char** argv) {
    gbbs::commandLine P(argc, argv, "Warm-start static k-median local search");

    auto app = [](auto& G, gbbs::commandLine P) {
        return gbbs::WarmStartStaticLocalSearch_runner(G, P);
    };

    return gbbs::kcenter_common::dispatch_main(argc, argv, app, /*mutates=*/false);
}
