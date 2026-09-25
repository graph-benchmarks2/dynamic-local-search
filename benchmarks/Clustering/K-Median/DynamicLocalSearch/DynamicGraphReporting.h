#pragma once

#include <iostream>

#include "benchmarks/Clustering/K-Median/DynamicLocalSearch/DynamicGraphDriverTypes.h"

namespace gbbs {
namespace kmedian {

inline void PrintDynamicGraphUpdate(
    const DynamicGraphUpdateResult& update_result) {

    std::cout << "### Dynamic update "
              << update_result.update_index
              << ": insert=("
              << update_result.u << ","
              << update_result.v << ")"
              << std::endl;

    std::cout << "num_distance_decreases = "
              << update_result.num_distance_decreases
              << std::endl;

    std::cout << "cost_before_update = "
              << update_result.cost_before_update
              << std::endl;

    std::cout << "cost_after_insertion_before_ls = "
              << update_result.cost_after_insertion_before_ls
              << std::endl;

    std::cout << "final_cost = "
              << update_result.final_cost
              << std::endl;

    std::cout << "triggered_swap = "
              << (update_result.triggered_swap ? "true" : "false")
              << std::endl;

    std::cout << "num_swaps = "
              << update_result.num_swaps
              << std::endl;

    std::cout << "num_static_recovery_swaps = "
              << update_result.num_static_swaps
              << std::endl;

    std::cout << "graph_rebuild_wall_time = "
              << update_result.graph_rebuild_wall_time
              << std::endl;

    std::cout << "distance_batch_wall_time = "
              << update_result.distance_batch_wall_time
              << std::endl;

    std::cout << "dynamic_wall_time = "
              << update_result.dynamic_wall_time
              << std::endl;

    std::cout << "distance_decrease_processing_wall_time = "
              << update_result.distance_decrease_processing_wall_time
              << std::endl;

    std::cout << "trigger_check_wall_time = "
              << update_result.trigger_check_wall_time
              << std::endl;

    std::cout << "first_swap_bookkeeping_wall_time = "
              << update_result.first_swap_bookkeeping_wall_time
              << std::endl;

    std::cout << "static_recovery_wall_time = "
              << update_result.static_recovery_wall_time
              << std::endl;

    std::cout << "dynamic_table_rebuild_wall_time = "
              << update_result.dynamic_table_rebuild_wall_time
              << std::endl;

    std::cout << "recovery_num_swap_table_builds = "
              << update_result.recovery_num_swap_table_builds
              << std::endl;

    std::cout << "recovery_initial_state_build_wall_time = "
              << update_result.recovery_initial_state_build_wall_time
              << std::endl;

    std::cout << "recovery_swap_table_build_wall_time = "
              << update_result.recovery_swap_table_build_wall_time
              << std::endl;

    std::cout << "recovery_state_rebuild_after_swap_wall_time = "
              << update_result.recovery_state_rebuild_after_swap_wall_time
              << std::endl;

    if (update_result.quality_evaluated) {
        std::cout << "cold_same_seed_cost = "
                << update_result.cold_same_seed_cost
                << std::endl;

        std::cout << "cold_best_cost = "
                << update_result.cold_best_cost
                << std::endl;

        std::cout << "maintained_over_cold_same_seed = "
                << update_result.maintained_over_cold_same_seed
                << std::endl;

        std::cout << "maintained_over_cold_best = "
                << update_result.maintained_over_cold_best
                << std::endl;

        std::cout << "cold_initializer_wall_time = "
                << update_result.cold_initializer_wall_time
                << std::endl;

        std::cout << "cold_local_search_wall_time = "
                << update_result.cold_local_search_wall_time
                << std::endl;

        std::cout << "cold_algorithm_wall_time = "
                << update_result.cold_initializer_wall_time +
                        update_result.cold_local_search_wall_time
                << std::endl;

        std::cout << "quality_baseline_wall_time = "
                << update_result.quality_baseline_wall_time
                << std::endl;
    }

    std::cout << "### ------------------------------------"
              << std::endl;
}

}  // namespace kmedian
}  // namespace gbbs
