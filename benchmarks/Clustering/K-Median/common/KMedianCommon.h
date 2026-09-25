#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <type_traits>
#include <utility>

#include "gbbs/gbbs.h"
#include "benchmarks/BFS/BFSdist/BFSdist.h"
#include "benchmarks/BFS/BFSdist/BFSdistSeq.h"
#include "benchmarks/PositiveWeightSSSP/DeltaStepping/DeltaStepping.h"
#include "benchmarks/PositiveWeightSSSP/Dijkstra/Dijkstra.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"

namespace gbbs {
    namespace kmedian {

        using Distance = float;

        constexpr uintE kNoCenter = UINT_E_MAX;

        inline Distance distance_infinity() {
            return std::numeric_limits<Distance>::max();
        }

        struct SSSPTiming {
            double total_time = 0.0;
            size_t calls = 0;

            void Add(double elapsed_seconds) {
                total_time += elapsed_seconds;
                ++calls;
            }
        };

        template <class RawDistance>
        inline bool is_unreachable(RawDistance d) {
            using R = typename std::decay<RawDistance>::type;

            if constexpr (std::is_floating_point<R>::value) {
                if (std::isinf(d)) {
                    return true;
                }

                if (d == static_cast<R>(std::numeric_limits<uintE>::max())) {
                    return true;
                }

                return d == std::numeric_limits<R>::max();
            } else {
                return d == std::numeric_limits<R>::max();
            }
        }

        struct KMedianSolutionState {
            parlay::sequence<uintE> centers;
            parlay::sequence<bool> is_center;

            parlay::sequence<uintE> nearest_center;
            parlay::sequence<uintE> second_nearest_center;

            parlay::sequence<Distance> nearest_dist;
            parlay::sequence<Distance> second_nearest_dist;

            double cost = 0.0;
            size_t unreachable_vertices = 0;
        };

        inline bool better_center(
            Distance candidate_dist,
            uintE candidate_center,
            Distance current_dist,
            uintE current_center) {

            if (candidate_dist < current_dist) {
                return true;
            }
            if (candidate_dist > current_dist) {
                return false;
            }

            return candidate_center < current_center;
        }

        template <class Graph>
        parlay::sequence<Distance> RunSSSP(
            Graph& G,
            uintE source,
            double delta,
            size_t num_buckets,
            bool single_core = false,
            bool sssp_verbose = false,
            SSSPTiming* sssp_timing = nullptr) {

            using W = typename Graph::weight_type;

            const size_t n = G.n;
            const Distance kInf = distance_infinity();

            timer sssp_timer;
            sssp_timer.start();

            parlay::sequence<Distance> result(n, kInf);

            if constexpr (std::is_same<W, gbbs::empty>::value) {
                auto distances = single_core
                    ? BFSdistSeq(G, source, std::numeric_limits<double>::infinity(), sssp_verbose)
                    : BFSdist(G, source, std::numeric_limits<double>::infinity(), sssp_verbose);

                parallel_for(0, n, [&](size_t v) {
                    auto raw = distances[v];

                    if (!is_unreachable(raw)) {
                        result[v] = static_cast<Distance>(raw);
                    }
                });
            } else {
                auto distances = single_core
                    ? Dijkstra(G, source, std::numeric_limits<double>::infinity(), sssp_verbose)
                    : DeltaStepping(
                        G,
                        source,
                        delta,
                        num_buckets,
                        std::numeric_limits<double>::infinity(),
                        sssp_verbose);

                parallel_for(0, n, [&](size_t v) {
                    auto raw = distances[v];

                    if (!is_unreachable(raw)) {
                        result[v] = static_cast<Distance>(raw);
                    }
                });
            }

            const double elapsed = sssp_timer.stop();

            if (sssp_timing != nullptr) {
                sssp_timing->Add(elapsed);
            }

            return result;
        }


        // Build the complete k-median solution state from an already available
        // exact distance metric. The metric must provide:
        //
        //     size_t n() const;
        //     distance(u, v)
        //
        // No graph traversal or SSSP is performed here.
        template <class Metric>
        KMedianSolutionState BuildSolutionStateFromDistances(
            const Metric& metric,
            const parlay::sequence<uintE>& centers) {

            const size_t n = metric.n();
            const Distance kInf = distance_infinity();

            KMedianSolutionState state;
            state.centers = centers;
            state.is_center = parlay::sequence<bool>(n, false);
            state.nearest_center = parlay::sequence<uintE>(n, kNoCenter);
            state.second_nearest_center = parlay::sequence<uintE>(n, kNoCenter);
            state.nearest_dist = parlay::sequence<Distance>(n, kInf);
            state.second_nearest_dist = parlay::sequence<Distance>(n, kInf);

            for (uintE c : centers) {
                if (c >= n) {
                    std::cout << "ERROR: Center " << c
                              << " is outside metric vertex range [0, " << n << ")."
                              << std::endl;
                    std::exit(-1);
                }

                if (state.is_center[c]) {
                    std::cout << "ERROR: Duplicate center " << c << "." << std::endl;
                    std::exit(-1);
                }

                state.is_center[c] = true;
            }

            for (uintE c : centers) {
                parallel_for(0, n, [&](size_t v) {
                    const auto raw = metric.distance(c, static_cast<uintE>(v));

                    if (is_unreachable(raw)) {
                        return;
                    }

                    const Distance d = static_cast<Distance>(raw);

                    if (better_center(
                            d,
                            c,
                            state.nearest_dist[v],
                            state.nearest_center[v])) {

                        state.second_nearest_dist[v] = state.nearest_dist[v];
                        state.second_nearest_center[v] = state.nearest_center[v];
                        state.nearest_dist[v] = d;
                        state.nearest_center[v] = c;
                    } else if (
                        c != state.nearest_center[v] &&
                        better_center(
                            d,
                            c,
                            state.second_nearest_dist[v],
                            state.second_nearest_center[v])) {

                        state.second_nearest_dist[v] = d;
                        state.second_nearest_center[v] = c;
                    }
                });
            }

            double cost = 0.0;
            size_t unreachable = 0;

            for (size_t v = 0; v < n; ++v) {
                if (state.nearest_center[v] == kNoCenter ||
                    state.nearest_dist[v] == kInf) {

                    ++unreachable;
                } else {
                    cost += static_cast<double>(state.nearest_dist[v]);
                }
            }

            state.unreachable_vertices = unreachable;
            state.cost =
                unreachable > 0
                    ? std::numeric_limits<double>::infinity()
                    : cost;

            return state;
        }

    }  // namespace kmedian
}  // namespace gbbs
