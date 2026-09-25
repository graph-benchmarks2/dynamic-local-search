#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>

#include "benchmarks/Clustering/K-Median/common/KMedianCommon.h"
#include "benchmarks/Clustering/K-Median/common/StaticSwapTable.h"

namespace gbbs {
    namespace kmedian {

        struct StaticLocalSearchResult {
            KMedianSolutionState state;
            StaticSwapTable final_table;

            size_t num_swaps = 0;
            size_t num_swap_table_builds = 0;

            double initial_cost = std::numeric_limits<double>::infinity();
            double final_cost = std::numeric_limits<double>::infinity();
            double final_max_gain = -std::numeric_limits<double>::infinity();
            double final_threshold = 0.0;

            double wall_time = 0.0;

            // Fine-grained local-search timing.
            double initial_state_build_wall_time = 0.0;
            double swap_table_build_wall_time = 0.0;
            double state_rebuild_after_swap_wall_time = 0.0;
        };

        template <class Metric>
        StaticLocalSearchResult RunStaticLocalSearchFromDistances(
            const Metric& metric,
            const parlay::sequence<uintE>& initial_centers,
            double ls_delta,
            bool verbose = false) {

            if (initial_centers.empty()) {
                std::cout << "ERROR: Static local search requires at least one center."
                          << std::endl;
                std::exit(-1);
            }

            if (ls_delta < 0.0) {
                std::cout << "ERROR: -ls-delta must be non-negative."
                          << std::endl;
                std::exit(-1);
            }

            const size_t k = initial_centers.size();

            timer wall_timer;
            wall_timer.start();

            timer initial_state_timer;
            initial_state_timer.start();

            auto state =
                BuildSolutionStateFromDistances(
                    metric,
                    initial_centers);

            StaticLocalSearchResult result;
            result.initial_state_build_wall_time =
                initial_state_timer.stop();
            result.initial_cost = state.cost;

            size_t iteration = 0;

            while (true) {
                timer table_timer;
                table_timer.start();

                StaticSwapTable table(
                    metric,
                    state);

                result.swap_table_build_wall_time +=
                    table_timer.stop();
                ++result.num_swap_table_builds;

                const double max_gain =
                    table.max_gain();

                const double threshold =
                    (ls_delta / static_cast<double>(k)) *
                    state.cost;

                if (verbose) {
                    std::cout << "### LS iteration " << iteration
                              << ": cost=" << state.cost
                              << " max_gain=" << max_gain
                              << " threshold=" << threshold
                              << std::endl;
                }

                result.final_max_gain = max_gain;
                result.final_threshold = threshold;

                if (!(max_gain > 0.0 && max_gain >= threshold)) {
                    result.final_table = std::move(table);
                    break;
                }

                const uintE x =
                    table.best_outgoing_center();

                const uintE y =
                    table.best_incoming_center();

                if (x == kNoCenter || y == kNoCenter) {
                    std::cout
                        << "ERROR: Swap table reports a significant gain "
                           "but no valid best swap."
                        << std::endl;
                    std::exit(-1);
                }

                if (verbose) {
                    std::cout << "### Performing swap "
                              << result.num_swaps + 1
                              << ": remove=" << x
                              << " add=" << y
                              << " gain=" << max_gain
                              << std::endl;
                }

                parlay::sequence<uintE> new_centers =
                    state.centers;

                bool found_x = false;

                for (size_t i = 0; i < new_centers.size(); ++i) {
                    if (new_centers[i] == x) {
                        new_centers[i] = y;
                        found_x = true;
                        break;
                    }
                }

                if (!found_x) {
                    std::cout << "ERROR: Best outgoing center "
                              << x
                              << " is not contained in the current center set."
                              << std::endl;
                    std::exit(-1);
                }

                const double old_cost =
                    state.cost;

                timer state_rebuild_timer;
                state_rebuild_timer.start();

                state =
                    BuildSolutionStateFromDistances(
                        metric,
                        new_centers);

                result.state_rebuild_after_swap_wall_time +=
                    state_rebuild_timer.stop();

                if (verbose) {
                    std::cout << "### Cost after swap: "
                              << state.cost
                              << " (decrease="
                              << old_cost - state.cost
                              << ")"
                              << std::endl;
                }

                ++result.num_swaps;
                ++iteration;
            }

            result.wall_time = wall_timer.stop();

            // Distances are already available. Therefore the entire measured
            // local search wall time is clustering time.

            result.final_cost = state.cost;
            result.state = std::move(state);

            return result;
        }

    }  // namespace kmedian
}  // namespace gbbs
