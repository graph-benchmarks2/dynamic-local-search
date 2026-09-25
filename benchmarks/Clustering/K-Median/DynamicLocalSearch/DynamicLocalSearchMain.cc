#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "gbbs/gbbs.h"

#include "benchmarks/Clustering/K-Center/common/DispatchMain.h"
#include "benchmarks/Clustering/K-Median/DynamicLocalSearch/DynamicGraphDriver.h"
#include "benchmarks/Clustering/K-Median/common/Initializers.h"

namespace gbbs {

parlay::sequence<uintE> ParseDynamicCenters(
    const std::string& text,
    size_t n) {

    std::vector<uintE> centers;

    if (text.empty()) {
        return parlay::sequence<uintE>();
    }

    std::stringstream stream(text);
    std::string token;

    while (std::getline(stream, token, ',')) {
        if (token.empty()) {
            std::cout << "ERROR: Empty entry in -centers."
                      << std::endl;
            std::exit(-1);
        }

        const uint64_t raw =
            static_cast<uint64_t>(
                std::stoull(token));

        if (raw >= n) {
            std::cout << "ERROR: Fixed center "
                      << raw
                      << " is outside [0," << n << ")."
                      << std::endl;
            std::exit(-1);
        }

        centers.push_back(
            static_cast<uintE>(raw));
    }

    return parlay::sequence<uintE>(
        centers.begin(),
        centers.end());
}

std::vector<size_t> ParseQualityCheckpoints(
    const std::string& text) {

    std::vector<size_t> checkpoints;

    if (text.empty()) {
        return checkpoints;
    }

    std::stringstream stream(text);
    std::string token;

    while (std::getline(stream, token, ',')) {
        if (token.empty()) {
            std::cout
                << "ERROR: Invalid -quality-cold-checkpoints specification: "
                << text
                << std::endl;
            std::exit(-1);
        }

        uint64_t raw = 0;

        try {
            raw = static_cast<uint64_t>(
                std::stoull(token));
        } catch (...) {
            std::cout
                << "ERROR: Invalid quality checkpoint: "
                << token
                << std::endl;
            std::exit(-1);
        }

        if (raw == 0) {
            std::cout
                << "ERROR: Quality checkpoints are 1-based and must be positive."
                << std::endl;
            std::exit(-1);
        }

        checkpoints.push_back(
            static_cast<size_t>(raw));
    }

    std::sort(
        checkpoints.begin(),
        checkpoints.end());

    checkpoints.erase(
        std::unique(
            checkpoints.begin(),
            checkpoints.end()),
        checkpoints.end());

    return checkpoints;
}

template <class Graph>
double DynamicLocalSearch_runner(
    Graph& G,
    commandLine P) {

    const size_t k =
        P.getOptionLongValue("-k", 10);

    const double delta =
        P.getOptionDoubleValue("-delta", 1.0);

    const double ls_delta =
        P.getOptionDoubleValue("-ls-delta", 0.1);

    const size_t num_buckets =
        P.getOptionLongValue("-nb", 32);

    const uint64_t seed =
        static_cast<uint64_t>(
            P.getOptionLongValue("-seed", 42));

    const std::string initializer =
        P.getOptionValue(
            "-init",
            "random");

    const std::string fixed_centers_text =
        P.getOptionValue(
            "-centers",
            "");

    const std::string updates_file =
        P.getOptionValue(
            "-updates",
            "");

    const bool single_core =
        P.getOption("-sc");

    const bool verbose =
        P.getOption("-verbose");

    const bool sssp_verbose =
        P.getOption("-sssp-verbose");

    const size_t cold_quality_trials =
        P.getOptionLongValue(
            "-quality-cold-trials",
            0);

    const uint64_t cold_quality_seed =
        static_cast<uint64_t>(
            P.getOptionLongValue(
                "-quality-cold-seed",
                static_cast<long>(seed)));

    const std::string cold_quality_checkpoints_text =
        P.getOptionValue(
            "-quality-cold-checkpoints",
            "");

    const std::vector<size_t> cold_quality_checkpoints =
        ParseQualityCheckpoints(
            cold_quality_checkpoints_text);

    if (updates_file.empty()) {
        std::cout << "ERROR: Dynamic driver requires -updates <file>."
                  << std::endl;
        std::exit(-1);
    }

    if (k == 0 || k > G.n) {
        std::cout << "ERROR: k must satisfy 1 <= k <= n."
                  << std::endl;
        std::exit(-1);
    }

    parlay::sequence<uintE> initial_centers;

    if (!fixed_centers_text.empty()) {
        initial_centers =
            ParseDynamicCenters(
                fixed_centers_text,
                G.n);

        if (initial_centers.size() != k) {
            std::cout << "ERROR: -centers contains "
                      << initial_centers.size()
                      << " centers but -k = " << k << "."
                      << std::endl;
            std::exit(-1);
        }
    } else {
        initial_centers =
            kmedian::InitializeCenters(
                G,
                initializer,
                k,
                delta,
                num_buckets,
                seed,
                single_core,
                sssp_verbose);
    }

    std::cout << "### Application: DynamicLocalSearch"
              << std::endl;

    std::cout << "### n: " << G.n
              << std::endl;

    std::cout << "### m: " << G.m
              << std::endl;

    std::cout << "### k: " << k
              << std::endl;

    std::cout << "### ls_delta: " << ls_delta
              << std::endl;

    std::cout << "### updates: " << updates_file
              << std::endl;

    if (cold_quality_trials > 0) {
        std::cout << "### quality_cold_trials: "
                  << cold_quality_trials
                  << std::endl;

        std::cout << "### quality_cold_seed: "
                  << cold_quality_seed
                  << std::endl;

        std::cout
            << "### quality baseline: fresh random cold starts"
            << std::endl;

        if (!cold_quality_checkpoints.empty()) {
            std::cout << "### quality_cold_checkpoints:";

            for (const size_t checkpoint :
                 cold_quality_checkpoints) {
                std::cout << " " << checkpoint;
            }

            std::cout << std::endl;
        }
    }

    std::cout << "### initial centers:";
    for (const uintE c : initial_centers) {
        std::cout << " " << c;
    }
    std::cout << std::endl;

    std::cout << "### ------------------------------------"
              << std::endl;

    auto result =
        kmedian::RunDynamicGraphSequence(
            G,
            std::move(initial_centers),
            updates_file,
            ls_delta,
            verbose,
            cold_quality_trials,
            cold_quality_seed,
            cold_quality_checkpoints);

    const size_t num_updates =
        result.updates.size();

    const double trigger_percentage =
        num_updates == 0
            ? 0.0
            : 100.0
                * static_cast<double>(
                    result.updates_triggering_swaps)
                / static_cast<double>(num_updates);

    std::cout << "### Dynamic summary"
              << std::endl;

    std::cout << "num_updates = "
              << num_updates
              << std::endl;

    std::cout << "total_distance_decreases = "
              << result.total_distance_decreases
              << std::endl;

    std::cout << "updates_triggering_swaps = "
              << result.updates_triggering_swaps
              << std::endl;

    std::cout << "updates_triggering_swaps_percent = "
              << trigger_percentage
              << std::endl;

    std::cout << "total_swaps = "
              << result.total_swaps
              << std::endl;

    std::cout << "final_cost = "
              << result.state.cost
              << std::endl;

    std::cout << "final_centers:";
    for (const uintE c : result.state.centers) {
        std::cout << " " << c;
    }
    std::cout << std::endl;

    std::cout << "initial_distance_matrix_wall_time = "
              << result.initial_distance_matrix_wall_time
              << std::endl;

    std::cout << "initial_distance_matrix_sssp_time = "
              << result.initial_distance_matrix_sssp_time
              << std::endl;

    std::cout << "initial_distance_matrix_sssp_calls = "
              << result.initial_distance_matrix_sssp_calls
              << std::endl;
              
    std::cout << "initial_local_search_wall_time = "
          << result.initial_local_search_wall_time
          << std::endl;

    std::cout << "initial_solution_setup_wall_time = "
            << result.initial_solution_setup_wall_time
            << std::endl;

    std::cout << "total_graph_rebuild_wall_time = "
              << result.total_graph_rebuild_wall_time
              << std::endl;

    std::cout << "total_distance_batch_wall_time = "
              << result.total_distance_batch_wall_time
              << std::endl;

    std::cout << "total_dynamic_wall_time = "
              << result.total_dynamic_wall_time
              << std::endl;

    std::cout << "total_distance_decrease_processing_wall_time = "
              << result.total_distance_decrease_processing_wall_time
              << std::endl;

    std::cout << "avg_distance_decrease_processing_time = "
              << (result.total_distance_decreases == 0
                      ? 0.0
                      : result.total_distance_decrease_processing_wall_time
                            / static_cast<double>(result.total_distance_decreases))
              << std::endl;

    std::cout << "total_trigger_check_wall_time = "
              << result.total_trigger_check_wall_time
              << std::endl;

    std::cout << "total_first_swap_bookkeeping_wall_time = "
              << result.total_first_swap_bookkeeping_wall_time
              << std::endl;

    std::cout << "total_static_recovery_wall_time = "
              << result.total_static_recovery_wall_time
              << std::endl;

    std::cout << "total_dynamic_table_rebuild_wall_time = "
              << result.total_dynamic_table_rebuild_wall_time
              << std::endl;

    std::cout << "total_recovery_num_swap_table_builds = "
              << result.total_recovery_num_swap_table_builds
              << std::endl;

    std::cout << "total_recovery_initial_state_build_wall_time = "
              << result.total_recovery_initial_state_build_wall_time
              << std::endl;

    std::cout << "total_recovery_swap_table_build_wall_time = "
              << result.total_recovery_swap_table_build_wall_time
              << std::endl;

    std::cout << "total_recovery_state_rebuild_after_swap_wall_time = "
              << result.total_recovery_state_rebuild_after_swap_wall_time
              << std::endl;

    // APSP maintenance and graph
    // rebuilding are timed separately and excluded from this benchmark time.
    return result.total_dynamic_wall_time;
}

}  // namespace gbbs

int main(int argc, char** argv) {
    auto app =
        [](auto& G, gbbs::commandLine P) {

            return gbbs::DynamicLocalSearch_runner(
                G,
                P);
        };

    return gbbs::kcenter_common::dispatch_main(
        argc,
        argv,
        app,
        /*mutates=*/false);
}
