#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "gbbs/gbbs.h"
#include "benchmarks/Clustering/K-Median/common/IncrementalAPSP.h"
#include "benchmarks/Clustering/K-Median/common/Initializers.h"
#include "benchmarks/Clustering/K-Median/StaticLocalSearch/StaticLocalSearch.h"
#include "parlay/sequence.h"

namespace gbbs {
    namespace kmedian {

        template <class W>
        struct EdgeInsertion {
            uintE u;
            uintE v;
            W w;
        };

        inline uint64_t UndirectedEdgeKey(uintE u, uintE v) {
            if (u > v) {
                std::swap(u, v);
            }

            return (static_cast<uint64_t>(u) << 32) | static_cast<uint64_t>(v);
        }

        template <class Graph>
        parlay::sequence<typename Graph::edge> ExtractUniqueUndirectedEdges(Graph& G) {
            using W = typename Graph::weight_type;
            using edge = typename Graph::edge;

            std::vector<edge> edges;
            edges.reserve(G.m / 2 + 1);

            for (uintE u = 0; u < G.n; ++u) {
                auto map_f = [&](const uintE& src, const uintE& v, const W& w) {
                    if (src < v) {
                        edges.emplace_back(src, v, w);
                    }
                };

                G.get_vertex(u).out_neighbors().map(map_f, false);
            }

            return parlay::sequence<edge>(edges.begin(), edges.end());
        }

        template <class Graph>
        parlay::sequence<EdgeInsertion<typename Graph::weight_type>> ReadInsertions(
            const std::string& filename,
            size_t n) {

            using W = typename Graph::weight_type;
            using insertion = EdgeInsertion<W>;

            std::ifstream in(filename);

            if (!in) {
                std::cout << "ERROR: Could not open insertion file: " << filename << std::endl;
                std::exit(-1);
            }

            std::vector<insertion> updates;

            while (true) {
                uint64_t u_raw;
                uint64_t v_raw;

                if (!(in >> u_raw >> v_raw)) {
                    break;
                }

                if (u_raw >= n || v_raw >= n) {
                    std::cout << "ERROR: Insertion endpoint outside graph: "
                              << u_raw << " " << v_raw
                              << std::endl;
                    std::exit(-1);
                }

                if (u_raw == v_raw) {
                    std::cout << "ERROR: Self-loop insertion is not supported: "
                              << u_raw << " " << v_raw
                              << std::endl;
                    std::exit(-1);
                }

                const uintE u = static_cast<uintE>(u_raw);
                const uintE v = static_cast<uintE>(v_raw);

                if constexpr (std::is_same<W, gbbs::empty>::value) {
                    updates.push_back(insertion{u, v, W{}});
                } else {
                    double w_raw;

                    if (!(in >> w_raw)) {
                        std::cout << "ERROR: Weighted insertion file must contain lines of the form 'u v w'."
                                  << std::endl;
                        std::exit(-1);
                    }

                    if (w_raw <= 0.0) {
                        std::cout << "ERROR: Edge weights must be positive." << std::endl;
                        std::exit(-1);
                    }

                    updates.push_back(insertion{u, v, static_cast<W>(w_raw)});
                }
            }

            return parlay::sequence<insertion>(updates.begin(), updates.end());
        }

        template <class W>
        using WarmTempGraph = gbbs::symmetric_graph<gbbs::symmetric_vertex, W>;

        template <class W>
        WarmTempGraph<W> BuildGraphFromEdges(
            const parlay::sequence<std::tuple<uintE, uintE, W>>& edges,
            size_t n) {

            return WarmTempGraph<W>::from_edges(edges, n);
        }

        struct WarmUpdateResult {
            size_t update_index = 0;

            uintE u = 0;
            uintE v = 0;

            double cost_before_update = 0.0;
            double cost_before_local_search = 0.0;
            double final_cost = 0.0;

            size_t num_swaps = 0;

            double update_wall_time = 0.0;

            // Local search timing inherited from the static
            // warm-start search on the updated metric.
            size_t num_swap_table_builds = 0;
            double initial_state_build_wall_time = 0.0;
            double swap_table_build_wall_time = 0.0;
            double state_rebuild_after_swap_wall_time = 0.0;

            // Optional quality cold-start baseline. These fields remain
            // untouched when cold_quality_trials == 0.
            double cold_same_seed_cost = 0.0;
            double cold_best_cost = 0.0;
            double maintained_over_cold_same_seed = 0.0;
            double maintained_over_cold_best = 0.0;

            // Timing for the cold baseline, kept separate. When multiple
            // cold trials are requested, these are sums across the trials at
            // this checkpoint.
            double cold_initializer_wall_time = 0.0;
            double cold_local_search_wall_time = 0.0;
            double quality_baseline_wall_time = 0.0;
            bool quality_evaluated = false;
        };

        template <class Graph>
        parlay::sequence<WarmUpdateResult> RunWarmStartSequence(
            Graph& initial_graph,
            IncrementalAPSP<Graph>& incremental_apsp,
            const parlay::sequence<uintE>& initial_centers,
            double initial_cost,
            const std::string& insertion_file,
            double ls_delta,
            bool verbose = false,
            size_t cold_quality_trials = 0,
            uint64_t cold_quality_seed = 42,
            const std::vector<size_t>& cold_quality_checkpoints = {}) {

            using W = typename Graph::weight_type;
            using edge = std::tuple<uintE, uintE, W>;

            const size_t n = initial_graph.n;

            auto base_edges = ExtractUniqueUndirectedEdges(initial_graph);
            parlay::sequence<edge> edges(base_edges.begin(), base_edges.end());

            std::unordered_set<uint64_t> present_edges;
            present_edges.reserve(edges.size() * 2 + 1);

            for (const auto& e : edges) {
                const uintE u = std::get<0>(e);
                const uintE v = std::get<1>(e);
                present_edges.insert(UndirectedEdgeKey(u, v));
            }

            auto updates = ReadInsertions<Graph>(insertion_file, n);
            parlay::sequence<uintE> centers = initial_centers;
            double current_cost = initial_cost;

            auto current_graph = BuildGraphFromEdges<W>(edges, n);

            parlay::sequence<WarmUpdateResult> results;
            results.reserve(updates.size());

            for (size_t i = 0; i < updates.size(); ++i) {
                const auto& upd = updates[i];
                const uintE u = upd.u;
                const uintE v = upd.v;

                const uint64_t key = UndirectedEdgeKey(u, v);
                if (present_edges.find(key) != present_edges.end()) {
                    std::cout << "ERROR: Update " << i
                              << " tries to insert an edge that already exists: "
                              << u << " " << v
                              << std::endl;
                    std::exit(-1);
                }

                auto decreases = incremental_apsp.GenerateBatch(
                    current_graph,
                    u,
                    v,
                    upd.w);

                for (const auto& decrease : decreases) {
                    incremental_apsp.set_distance(
                        decrease.p,
                        decrease.q,
                        decrease.new_distance);
                }

                auto ls_result = RunStaticLocalSearchFromDistances(
                    incremental_apsp,
                    centers,
                    ls_delta,
                    verbose);

                WarmUpdateResult update_result;
                update_result.update_index = i;
                update_result.u = u;
                update_result.v = v;
                update_result.cost_before_update = current_cost;
                update_result.cost_before_local_search = ls_result.initial_cost;
                update_result.final_cost = ls_result.final_cost;
                update_result.num_swaps = ls_result.num_swaps;
                update_result.update_wall_time = ls_result.wall_time;
                update_result.num_swap_table_builds =
                    ls_result.num_swap_table_builds;
                update_result.initial_state_build_wall_time =
                    ls_result.initial_state_build_wall_time;
                update_result.swap_table_build_wall_time =
                    ls_result.swap_table_build_wall_time;
                update_result.state_rebuild_after_swap_wall_time =
                    ls_result.state_rebuild_after_swap_wall_time;


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

                    double best_cold_cost = std::numeric_limits<double>::infinity();
                    double same_seed_cost = std::numeric_limits<double>::infinity();

                    for (size_t trial = 0; trial < cold_quality_trials; ++trial) {
                        const uint64_t trial_seed =
                            cold_quality_seed + static_cast<uint64_t>(trial);

                        timer cold_initializer_timer;
                        cold_initializer_timer.start();

                        const auto cold_centers =
                            RandomCenters(n, centers.size(), trial_seed);

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
                            same_seed_cost = cold_result.final_cost;
                        }

                        best_cold_cost =
                            std::min(best_cold_cost, cold_result.final_cost);
                    }

                    update_result.quality_baseline_wall_time =
                        quality_timer.stop();
                    update_result.cold_same_seed_cost = same_seed_cost;
                    update_result.cold_best_cost = best_cold_cost;
                    update_result.maintained_over_cold_same_seed =
                        update_result.final_cost / same_seed_cost;
                    update_result.maintained_over_cold_best =
                        update_result.final_cost / best_cold_cost;
                }

                centers = ls_result.state.centers;
                current_cost = ls_result.final_cost;
                results.push_back(update_result);

                present_edges.insert(key);
                edges.push_back(std::make_tuple(u, v, upd.w));
                current_graph = BuildGraphFromEdges<W>(edges, n);

                std::cout << "### Warm update " << i << ": insert=(" << u << "," << v << ")" << std::endl;
                std::cout << "num_distance_decreases = " << decreases.size() << std::endl;
                std::cout << "cost_before_update = " << update_result.cost_before_update << std::endl;
                std::cout << "cost_after_insertion_before_ls = "
                          << update_result.cost_before_local_search
                          << std::endl;
                std::cout << "final_cost = " << update_result.final_cost << std::endl;
                std::cout << "num_swaps = " << update_result.num_swaps << std::endl;
                std::cout << "update_wall_time = "
                          << update_result.update_wall_time
                          << std::endl;
                std::cout << "num_swap_table_builds = "
                          << update_result.num_swap_table_builds
                          << std::endl;
                std::cout << "initial_state_build_wall_time = "
                          << update_result.initial_state_build_wall_time
                          << std::endl;
                std::cout << "swap_table_build_wall_time = "
                          << update_result.swap_table_build_wall_time
                          << std::endl;
                std::cout << "state_rebuild_after_swap_wall_time = "
                          << update_result.state_rebuild_after_swap_wall_time
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

                std::cout << "### ------------------------------------" << std::endl;
            }

            return results;
        }

    }  // namespace kmedian
}  // namespace gbbs
