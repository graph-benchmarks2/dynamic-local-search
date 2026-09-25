#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "gbbs/gbbs.h"

#include "benchmarks/Clustering/K-Median/DynamicLocalSearch/DynamicLocalSearch.h"
#include "benchmarks/Clustering/K-Median/common/IncrementalAPSP.h"
#include "benchmarks/Clustering/K-Median/DynamicLocalSearch/DynamicGraphDriverTypes.h"
#include "benchmarks/Clustering/K-Median/DynamicLocalSearch/DynamicGraphReporting.h"
#include "benchmarks/Clustering/K-Median/StaticLocalSearch/StaticLocalSearch.h"
#include "benchmarks/Clustering/K-Median/DynamicLocalSearch/DynamicGraphIO.h"
#include "benchmarks/Clustering/K-Median/common/StaticSwapTable.h"
#include "benchmarks/Clustering/K-Median/common/Initializers.h"
#include "parlay/sequence.h"

namespace gbbs {
namespace kmedian {

template <class W>
double DynamicInsertionWeight(const W& weight) {
    if constexpr (std::is_same<W, gbbs::empty>::value) {
        return 1.0;
    } else {
        return static_cast<double>(weight);
    }
}

// Correctness-oriented GBBS driver.
//   * initial all-pairs distances are obtained by repeated GBBS SSSP;
//   * one insertion is converted into a batch of exact distance decreases;
template <class Graph>
DynamicGraphSequenceResult RunDynamicGraphSequence(
    Graph& initial_graph,
    parlay::sequence<uintE> initial_centers,
    const std::string& insertion_file,
    double ls_delta,
    bool verbose = false,
    size_t cold_quality_trials = 0,
    uint64_t cold_quality_seed = 42,
    const std::vector<size_t>& cold_quality_checkpoints = {}) {

    using W = typename Graph::weight_type;
    using edge = std::tuple<uintE, uintE, W>;

    if (initial_centers.empty()) {
        std::cout << "ERROR: Dynamic graph driver requires a nonempty center set."
                  << std::endl;
        std::exit(-1);
    }

    const size_t n = initial_graph.n;

    auto base_edges = ExtractDynamicUniqueUndirectedEdges(initial_graph);
    parlay::sequence<edge> edges(base_edges.begin(), base_edges.end());

    std::unordered_set<uint64_t> present_edges;
    present_edges.reserve(edges.size() * 2 + 1);

    for (const auto& e : edges) {
        present_edges.insert(
            DynamicUndirectedEdgeKey(
                std::get<0>(e),
                std::get<1>(e)));
    }

    const auto insertions =
        ReadDynamicInsertions<Graph>(
            insertion_file,
            n);

    // Use one canonical simple symmetric graph type for incremental APSP
    // maintenance. 
    auto current_graph =
        BuildDynamicGraphFromEdges<W>(
            edges,
            n);

    DynamicGraphSequenceResult sequence_result;
    sequence_result.updates.reserve(insertions.size());

    SSSPTiming initial_distance_sssp;

    // Full initial-solution setup: exact APSP + initial static local search.
    timer initial_solution_setup_timer;
    initial_solution_setup_timer.start();

    timer initial_distance_timer;
    initial_distance_timer.start();

    IncrementalAPSP<decltype(current_graph)> incremental_apsp(
        current_graph,
        std::numeric_limits<double>::infinity(),
        &initial_distance_sssp);

    sequence_result.initial_distance_matrix_wall_time =
        initial_distance_timer.stop();

    sequence_result.initial_distance_matrix_sssp_time =
        initial_distance_sssp.total_time;

    sequence_result.initial_distance_matrix_sssp_calls =
        initial_distance_sssp.calls;

    timer initial_local_search_timer;
    initial_local_search_timer.start();

    auto initial_ls =
        RunStaticLocalSearchFromDistances(
            incremental_apsp,
            initial_centers,
            ls_delta,
            verbose);

    sequence_result.initial_local_search_wall_time =
        initial_local_search_timer.stop();

    sequence_result.initial_solution_setup_wall_time =
        initial_solution_setup_timer.stop();

    sequence_result.state =
        std::move(initial_ls.state);

    std::cout << "### Initial locally stable solution"
              << std::endl;

    std::cout << "initial_final_cost = "
              << initial_ls.final_cost
              << std::endl;

    std::cout << "initial_num_swaps = "
              << initial_ls.num_swaps
              << std::endl;

    std::cout << "initial_centers:";
    for (const uintE c : sequence_result.state.centers) {
        std::cout << " " << c;
    }
    std::cout << std::endl;

    std::cout << "### ===================================="
              << std::endl;

    DynamicSwapTable sw_decr(
        initial_ls.final_table,
        n);


    for (size_t i = 0; i < insertions.size(); ++i) {
        const auto& insertion = insertions[i];

        const uintE u = insertion.u;
        const uintE v = insertion.v;
        const uint64_t key = DynamicUndirectedEdgeKey(u, v);

        if (present_edges.find(key) != present_edges.end()) {
            std::cout << "ERROR: Dynamic update " << i
                      << " tries to insert an existing edge: "
                      << u << " " << v
                      << std::endl;
            std::exit(-1);
        }

        DynamicGraphUpdateResult update_result;
        update_result.update_index = i;
        update_result.u = u;
        update_result.v = v;
        update_result.cost_before_update =
            sequence_result.state.cost;


        timer distance_batch_timer;
        distance_batch_timer.start();

        auto batch = incremental_apsp.GenerateBatch(
            current_graph,
            u,
            v,
            insertion.w);

        update_result.distance_batch_wall_time =
            distance_batch_timer.stop();

        update_result.num_distance_decreases =
            batch.size();


        present_edges.insert(key);
        edges.push_back(
            std::make_tuple(
                u,
                v,
                insertion.w));

        timer graph_rebuild_timer;
        graph_rebuild_timer.start();

        auto updated_graph =
            BuildDynamicGraphFromEdges<W>(
                edges,
                n);

        update_result.graph_rebuild_wall_time =
            graph_rebuild_timer.stop();


        auto get_distance =
            [&](uintE a, uintE b) {

                return incremental_apsp.distance(a, b);
            };

        auto set_distance =
            [&](uintE a, uintE b, double value) {

                incremental_apsp.set_distance(
                    a,
                    b,
                    value);
            };

        auto static_recovery =
            [&](const parlay::sequence<uintE>& swapped_centers) {

                return RunStaticLocalSearchFromDistances(
                    incremental_apsp,
                    swapped_centers,
                    ls_delta,
                    verbose);
            };

        timer dynamic_timer;
        dynamic_timer.start();

        const auto batch_result =
            ProcessDistanceDecreaseBatch(
                sw_decr,
                sequence_result.state,
                batch,
                ls_delta,
                get_distance,
                set_distance,
                static_recovery,
                verbose);

        update_result.dynamic_wall_time =
            dynamic_timer.stop();

        update_result.cost_after_insertion_before_ls =
            batch_result.cost_after_batch_before_ls;

        update_result.final_cost =
            batch_result.final_cost;

        update_result.num_swaps =
            batch_result.num_swaps;

        update_result.num_static_swaps =
            batch_result.num_static_swaps;

        update_result.triggered_swap =
            batch_result.triggered_swap;

        update_result.distance_decrease_processing_wall_time =
            batch_result.distance_decrease_processing_wall_time;
        update_result.trigger_check_wall_time =
            batch_result.trigger_check_wall_time;
        update_result.first_swap_bookkeeping_wall_time =
            batch_result.first_swap_bookkeeping_wall_time;
        update_result.static_recovery_wall_time =
            batch_result.static_recovery_wall_time;
        update_result.dynamic_table_rebuild_wall_time =
            batch_result.dynamic_table_rebuild_wall_time;

        update_result.recovery_num_swap_table_builds =
            batch_result.recovery_num_swap_table_builds;
        update_result.recovery_initial_state_build_wall_time =
            batch_result.recovery_initial_state_build_wall_time;
        update_result.recovery_swap_table_build_wall_time =
            batch_result.recovery_swap_table_build_wall_time;
        update_result.recovery_state_rebuild_after_swap_wall_time =
            batch_result.recovery_state_rebuild_after_swap_wall_time;


        const bool evaluate_quality =
            cold_quality_trials > 0 &&
            (cold_quality_checkpoints.empty() ||
            std::binary_search(
                cold_quality_checkpoints.begin(),
                cold_quality_checkpoints.end(),
                i + 1));

        if (evaluate_quality) {
            update_result.quality_evaluated = true;

            timer quality_timer;
            quality_timer.start();

            double best_cold_cost =
                std::numeric_limits<double>::infinity();

            double same_seed_cost =
                std::numeric_limits<double>::infinity();

            for (size_t trial = 0;
                trial < cold_quality_trials;
                ++trial) {

                const uint64_t trial_seed =
                    cold_quality_seed +
                    static_cast<uint64_t>(trial);

                timer cold_initializer_timer;
                cold_initializer_timer.start();

                const auto cold_centers =
                    RandomCenters(
                        n,
                        sequence_result.state.centers.size(),
                        trial_seed);

                update_result.cold_initializer_wall_time +=
                    cold_initializer_timer.stop();

                const auto cold_result =
                    RunStaticLocalSearchFromDistances(
                        incremental_apsp,
                        cold_centers,
                        ls_delta,
                        /*verbose=*/false);

                update_result.cold_local_search_wall_time +=
                    cold_result.wall_time;

                if (trial == 0) {
                    same_seed_cost =
                        cold_result.final_cost;
                }

                best_cold_cost =
                    std::min(
                        best_cold_cost,
                        cold_result.final_cost);
            }

            update_result.quality_baseline_wall_time =
                quality_timer.stop();

            update_result.cold_same_seed_cost =
                same_seed_cost;

            update_result.cold_best_cost =
                best_cold_cost;

            update_result.maintained_over_cold_same_seed =
                update_result.final_cost /
                same_seed_cost;

            update_result.maintained_over_cold_best =
                update_result.final_cost /
                best_cold_cost;
        }

        sequence_result.total_distance_decreases +=
            update_result.num_distance_decreases;

        sequence_result.total_swaps +=
            update_result.num_swaps;

        if (update_result.triggered_swap) {
            ++sequence_result.updates_triggering_swaps;
        }

        sequence_result.total_graph_rebuild_wall_time +=
            update_result.graph_rebuild_wall_time;

        sequence_result.total_distance_batch_wall_time +=
            update_result.distance_batch_wall_time;

        sequence_result.total_dynamic_wall_time +=
            update_result.dynamic_wall_time;

        sequence_result.total_distance_decrease_processing_wall_time +=
            update_result.distance_decrease_processing_wall_time;
        sequence_result.total_trigger_check_wall_time +=
            update_result.trigger_check_wall_time;
        sequence_result.total_first_swap_bookkeeping_wall_time +=
            update_result.first_swap_bookkeeping_wall_time;
        sequence_result.total_static_recovery_wall_time +=
            update_result.static_recovery_wall_time;
        sequence_result.total_dynamic_table_rebuild_wall_time +=
            update_result.dynamic_table_rebuild_wall_time;

        sequence_result.total_recovery_num_swap_table_builds +=
            update_result.recovery_num_swap_table_builds;
        sequence_result.total_recovery_initial_state_build_wall_time +=
            update_result.recovery_initial_state_build_wall_time;
        sequence_result.total_recovery_swap_table_build_wall_time +=
            update_result.recovery_swap_table_build_wall_time;
        sequence_result.total_recovery_state_rebuild_after_swap_wall_time +=
            update_result.recovery_state_rebuild_after_swap_wall_time;

        sequence_result.updates.push_back(
            update_result);

        PrintDynamicGraphUpdate(update_result);

        // The next APSP repair must see all insertions processed so far.
        // current_graph and updated_graph have the same canonical graph type.
        current_graph = std::move(updated_graph);
    }

    return sequence_result;
}

}  // namespace kmedian
}  // namespace gbbs
