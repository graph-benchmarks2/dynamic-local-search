#pragma once

#include <cstddef>

#include "benchmarks/Clustering/K-Median/common/KMedianCommon.h"
#include "parlay/sequence.h"

namespace gbbs {
namespace kmedian {

struct DynamicGraphUpdateResult {
    size_t update_index = 0;
    uintE u = 0;
    uintE v = 0;

    size_t num_distance_decreases = 0;
    size_t num_swaps = 0;
    size_t num_static_swaps = 0;

    bool triggered_swap = false;

    double cost_before_update = 0.0;
    double cost_after_insertion_before_ls = 0.0;
    double final_cost = 0.0;

    double graph_rebuild_wall_time = 0.0;
    double distance_batch_wall_time = 0.0;

    double dynamic_wall_time = 0.0;

    // Fine-grained dynamic local-search timing.
    double distance_decrease_processing_wall_time = 0.0;
    double trigger_check_wall_time = 0.0;
    double first_swap_bookkeeping_wall_time = 0.0;
    double static_recovery_wall_time = 0.0;
    double dynamic_table_rebuild_wall_time = 0.0;

    size_t recovery_num_swap_table_builds = 0;
    double recovery_initial_state_build_wall_time = 0.0;
    double recovery_swap_table_build_wall_time = 0.0;
    double recovery_state_rebuild_after_swap_wall_time = 0.0;

    // Optional cold-start quality baseline.
    double cold_same_seed_cost = 0.0;
    double cold_best_cost = 0.0;

    double maintained_over_cold_same_seed = 0.0;
    double maintained_over_cold_best = 0.0;

    double cold_initializer_wall_time = 0.0;
    double cold_local_search_wall_time = 0.0;
    double quality_baseline_wall_time = 0.0;

    bool quality_evaluated = false;
};

struct DynamicGraphSequenceResult {
    KMedianSolutionState state;
    parlay::sequence<DynamicGraphUpdateResult> updates;

    size_t total_distance_decreases = 0;
    size_t total_swaps = 0;
    size_t updates_triggering_swaps = 0;

    double initial_distance_matrix_wall_time = 0.0;
    double initial_distance_matrix_sssp_time = 0.0;
    size_t initial_distance_matrix_sssp_calls = 0;

    double initial_local_search_wall_time = 0.0;
    double initial_solution_setup_wall_time = 0.0;

    double total_graph_rebuild_wall_time = 0.0;
    double total_distance_batch_wall_time = 0.0;

    double total_dynamic_wall_time = 0.0;

    // Aggregated fine-grained dynamic local-search timing.
    double total_distance_decrease_processing_wall_time = 0.0;
    double total_trigger_check_wall_time = 0.0;
    double total_first_swap_bookkeeping_wall_time = 0.0;
    double total_static_recovery_wall_time = 0.0;
    double total_dynamic_table_rebuild_wall_time = 0.0;

    size_t total_recovery_num_swap_table_builds = 0;
    double total_recovery_initial_state_build_wall_time = 0.0;
    double total_recovery_swap_table_build_wall_time = 0.0;
    double total_recovery_state_rebuild_after_swap_wall_time = 0.0;

};


}  // namespace kmedian
}  // namespace gbbs
