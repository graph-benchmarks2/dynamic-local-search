#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>

#include "benchmarks/Clustering/K-Median/common/DynamicSwapTable.h"
#include "benchmarks/Clustering/K-Median/common/KMedianCommon.h"
#include "benchmarks/Clustering/K-Median/common/DistanceDecrease.h"
#include "parlay/sequence.h"

namespace gbbs {
namespace kmedian {

struct DynamicBatchResult {
    size_t num_distance_decreases = 0;
    size_t num_swaps = 0;
    size_t num_static_swaps = 0;

    bool triggered_swap = false;

    uintE first_outgoing_center = kNoCenter;
    uintE first_incoming_center = kNoCenter;

    double cost_before_batch = std::numeric_limits<double>::infinity();
    double cost_after_batch_before_ls = std::numeric_limits<double>::infinity();
    double cost_after_first_swap = std::numeric_limits<double>::infinity();
    double final_cost = std::numeric_limits<double>::infinity();

    double trigger_max_gain = -std::numeric_limits<double>::infinity();
    double trigger_threshold = 0.0;

    double final_max_gain = -std::numeric_limits<double>::infinity();
    double final_threshold = 0.0;

    // Fine-grained timing for the dynamic local-search work.
    double distance_decrease_processing_wall_time = 0.0;
    double trigger_check_wall_time = 0.0;
    double first_swap_bookkeeping_wall_time = 0.0;
    double static_recovery_wall_time = 0.0;
    double dynamic_table_rebuild_wall_time = 0.0;

    // Fine-grained breakdown of the static recovery, when one is triggered.
    size_t recovery_num_swap_table_builds = 0;
    double recovery_initial_state_build_wall_time = 0.0;
    double recovery_swap_table_build_wall_time = 0.0;
    double recovery_state_rebuild_after_swap_wall_time = 0.0;
};

inline parlay::sequence<uintE> CentersAfterSwap(
    const parlay::sequence<uintE>& centers,
    uintE outgoing_center,
    uintE incoming_center) {

    parlay::sequence<uintE> swapped = centers;
    bool found = false;

    for (size_t i = 0; i < swapped.size(); ++i) {
        if (swapped[i] == outgoing_center) {
            swapped[i] = incoming_center;
            found = true;
            break;
        }
    }

    if (!found) {
        std::cout << "ERROR: Dynamic local search could not find the outgoing center."
                  << std::endl;
        std::exit(-1);
    }

    return swapped;
}

// Processes one adversarial batch of single-pair distance decreases.
template <
    class GetDistance,
    class SetDistance,
    class StaticRecovery>
DynamicBatchResult ProcessDistanceDecreaseBatch(
    DynamicSwapTable& sw_decr,
    KMedianSolutionState& state,
    const parlay::sequence<DistanceDecrease>& batch,
    double ls_delta,
    GetDistance get_distance,
    SetDistance set_distance,
    StaticRecovery static_recovery,
    bool verbose = false) {

    if (state.centers.empty()) {
        std::cout << "ERROR: Dynamic local search requires at least one center."
                  << std::endl;
        std::exit(-1);
    }

    if (ls_delta < 0.0) {
        std::cout << "ERROR: ls_delta must be non-negative."
                  << std::endl;
        std::exit(-1);
    }

    const size_t n = state.is_center.size();
    const size_t k = state.centers.size();

    DynamicBatchResult result;
    result.num_distance_decreases = batch.size();
    result.cost_before_batch = state.cost;

    timer decrease_processing_timer;
    decrease_processing_timer.start();

    for (const auto& decrease : batch) {
        const uintE p = decrease.p;
        const uintE q = decrease.q;

        if (p == q || p >= n || q >= n) {
            std::cout << "ERROR: Invalid distance decrease in dynamic batch."
                      << std::endl;
            std::exit(-1);
        }

        const double old_value = static_cast<double>(get_distance(p, q));
        const double new_value = decrease.new_distance;

        if (!(new_value < old_value)) {
            std::cout << "ERROR: Dynamic batch contains a non-decreasing update: "
                      << "(" << p << "," << q << ") "
                      << old_value << " -> " << new_value
                      << std::endl;
            std::exit(-1);
        }

        auto old_distance = [&](uintE u, uintE v) {
            return get_distance(u, v);
        };

        auto new_distance = [&](uintE u, uintE v) {
            if ((u == p && v == q) || (u == q && v == p)) {
                return new_value;
            }
            return static_cast<double>(get_distance(u, v));
        };

        sw_decr.ApplyDistanceDecrease(
            p, q, state, old_distance, new_distance);

        set_distance(p, q, new_value);
    }

    result.distance_decrease_processing_wall_time =
        decrease_processing_timer.stop();

    result.cost_after_batch_before_ls = state.cost;

    timer trigger_check_timer;
    trigger_check_timer.start();
    result.trigger_max_gain = sw_decr.max_gain();
    result.trigger_threshold =
        (ls_delta / static_cast<double>(k)) * state.cost;
    result.trigger_check_wall_time =
        trigger_check_timer.stop();

    if (verbose) {
        std::cout << "### Dynamic batch: cost=" << state.cost
                  << " max_gain=" << result.trigger_max_gain
                  << " threshold=" << result.trigger_threshold
                  << std::endl;
    }

    if (!(result.trigger_max_gain > 0.0
          && result.trigger_max_gain >= result.trigger_threshold)) {

        result.final_cost = state.cost;
        result.final_max_gain = result.trigger_max_gain;
        result.final_threshold = result.trigger_threshold;
        return result;
    }

    result.triggered_swap = true;
    result.first_outgoing_center = sw_decr.best_outgoing_center();
    result.first_incoming_center = sw_decr.best_incoming_center();

    if (result.first_outgoing_center == kNoCenter
        || result.first_incoming_center == kNoCenter) {

        std::cout << "ERROR: SW_decr reports a significant gain but no valid swap."
                  << std::endl;
        std::exit(-1);
    }

    if (verbose) {
        std::cout << "### Dynamic first swap: remove="
                  << result.first_outgoing_center
                  << " add=" << result.first_incoming_center
                  << " gain=" << result.trigger_max_gain
                  << std::endl;
    }

    timer first_swap_timer;
    first_swap_timer.start();

    const auto swapped_centers = CentersAfterSwap(
        state.centers,
        result.first_outgoing_center,
        result.first_incoming_center);

    result.first_swap_bookkeeping_wall_time =
        first_swap_timer.stop();

    // SW_decr is stale from this point onward.
    timer recovery_timer;
    recovery_timer.start();

    auto recovery = static_recovery(swapped_centers);

    result.static_recovery_wall_time =
        recovery_timer.stop();
    result.recovery_num_swap_table_builds =
        recovery.num_swap_table_builds;
    result.recovery_initial_state_build_wall_time =
        recovery.initial_state_build_wall_time;
    result.recovery_swap_table_build_wall_time =
        recovery.swap_table_build_wall_time;
    result.recovery_state_rebuild_after_swap_wall_time =
        recovery.state_rebuild_after_swap_wall_time;

    result.cost_after_first_swap = recovery.initial_cost;
    result.num_static_swaps = recovery.num_swaps;
    result.num_swaps = 1 + recovery.num_swaps;
    result.final_cost = recovery.final_cost;

    const double expected_first_swap_cost =
        result.cost_after_batch_before_ls - result.trigger_max_gain;

    const double scale = std::max(
        1.0,
        std::max(
            std::abs(expected_first_swap_cost),
            std::abs(result.cost_after_first_swap)));

    if (std::abs(
            expected_first_swap_cost - result.cost_after_first_swap)
        > 1e-6 * scale) {

        std::cout << "ERROR: First dynamic swap gain disagrees with the "
                     "recomputed post-swap cost."
                  << std::endl;
        std::exit(-1);
    }

    state = std::move(recovery.state);

    // Static recovery terminated only after constructing the locally stable
    // SW_stat for its final state.
    timer table_rebuild_timer;
    table_rebuild_timer.start();

    sw_decr.Build(recovery.final_table, n);

    result.dynamic_table_rebuild_wall_time =
        table_rebuild_timer.stop();

    result.final_max_gain = sw_decr.max_gain();
    result.final_threshold =
        (ls_delta / static_cast<double>(k)) * state.cost;

    if (verbose) {
        std::cout << "### Dynamic batch result: final_cost="
                  << result.final_cost
                  << " total_swaps=" << result.num_swaps
                  << " final_max_gain=" << result.final_max_gain
                  << " final_threshold=" << result.final_threshold
                  << std::endl;
    }

    return result;
}

}  // namespace kmedian
}  // namespace gbbs
