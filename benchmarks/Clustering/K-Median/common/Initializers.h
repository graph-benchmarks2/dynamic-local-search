#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "benchmarks/Clustering/K-Center/Gonzalez/Gonzalez.h"
#include "benchmarks/Clustering/K-Median/common/KMedianCommon.h"
#include "parlay/sequence.h"

namespace gbbs {
    namespace kmedian {

        inline parlay::sequence<uintE> RandomCenters(size_t n, size_t k, uint64_t seed) {
            if (k > n) {
                std::cout << "ERROR: k = " << k
                          << " exceeds number of vertices n = " << n << "."
                          << std::endl;
                std::exit(-1);
            }

            std::vector<uintE> vertices(n);
            std::iota(vertices.begin(), vertices.end(), uintE{0});

            std::mt19937_64 rng(seed);

            for (size_t i = 0; i < k; ++i) {
                std::uniform_int_distribution<size_t> pick(i, n - 1);
                size_t j = pick(rng);
                std::swap(vertices[i], vertices[j]);
            }

            parlay::sequence<uintE> centers(k);
            for (size_t i = 0; i < k; ++i) {
                centers[i] = vertices[i];
            }

            return centers;
        }

        template <class Graph>
        parlay::sequence<uintE> GonzalezCenters(
            Graph& G,
            size_t k,
            double delta,
            size_t num_buckets,
            uint64_t seed,
            bool single_core = false,
            bool sssp_verbose = false,
            SSSPTiming* sssp_timing = nullptr) {

            if (k > G.n) {
                std::cout << "ERROR: k = " << k
                          << " exceeds number of vertices n = " << G.n << "."
                          << std::endl;
                std::exit(-1);
            }

            timer distance_timer;
            if (sssp_timing != nullptr) {
                distance_timer.start();
            }

            auto result = gbbs::Gonzalez(
                G, static_cast<uintE>(k), delta, num_buckets, seed, single_core, sssp_verbose);

            if (sssp_timing != nullptr) {
                sssp_timing->total_time += distance_timer.stop();
            }

            return result.centers;
        }

        template <class Graph>
        parlay::sequence<uintE> GreedyKMedianCenters(
            Graph& G,
            size_t k,
            double delta,
            size_t num_buckets,
            uint64_t seed,
            bool single_core = false,
            bool sssp_verbose = false,
            SSSPTiming* sssp_timing = nullptr) {

            const size_t n = G.n;
            const Distance kInf = distance_infinity();

            if (k == 0) {
                return {};
            }

            if (k > n) {
                std::cout << "ERROR: k = " << k
                          << " exceeds number of vertices n = " << n << "."
                          << std::endl;
                std::exit(-1);
            }

            const uintE first = RandomCenters(n, 1, seed)[0];

            parlay::sequence<uintE> centers;
            centers.push_back(first);

            parlay::sequence<bool> is_center(n, false);
            is_center[first] = true;

            auto dist_to_centers = RunSSSP(
                G, first, delta, num_buckets, single_core, sssp_verbose, sssp_timing);

            for (size_t round = 1; round < k; ++round) {
                double best_cost = std::numeric_limits<double>::infinity();
                uintE best_candidate = kNoCenter;
                parlay::sequence<Distance> best_distances;

                for (uintE y = 0; y < n; ++y) {
                    if (is_center[y]) {
                        continue;
                    }

                    auto dist_y = RunSSSP(
                        G, y, delta, num_buckets, single_core, sssp_verbose, sssp_timing);

                    double candidate_cost = 0.0;
                    bool unreachable = false;

                    for (size_t v = 0; v < n; ++v) {
                        const Distance d = std::min(dist_to_centers[v], dist_y[v]);
                        if (d == kInf) {
                            unreachable = true;
                            break;
                        }

                        candidate_cost += static_cast<double>(d);
                    }

                    if (unreachable) {
                        candidate_cost = std::numeric_limits<double>::infinity();
                    }

                    if (candidate_cost < best_cost ||
                        (candidate_cost == best_cost && y < best_candidate)) {

                        best_cost = candidate_cost;
                        best_candidate = y;
                        best_distances = std::move(dist_y);
                    }
                }

                if (best_candidate == kNoCenter) {
                    std::cout << "ERROR: Greedy k-median initializer could not find another center."
                              << std::endl;
                    std::exit(-1);
                }

                centers.push_back(best_candidate);
                is_center[best_candidate] = true;

                parallel_for(0, n, [&](size_t v) {
                    dist_to_centers[v] = std::min(dist_to_centers[v], best_distances[v]);
                });
            }

            return centers;
        }

        template <class Graph>
        parlay::sequence<uintE> InitializeCenters(
            Graph& G,
            const std::string& mode,
            size_t k,
            double delta,
            size_t num_buckets,
            uint64_t seed,
            bool single_core = false,
            bool sssp_verbose = false,
            SSSPTiming* sssp_timing = nullptr) {

            if (mode == "greedy") {
                return GreedyKMedianCenters(
                    G, k, delta, num_buckets, seed, single_core, sssp_verbose, sssp_timing);
            }

            if (mode == "random") {
                return RandomCenters(G.n, k, seed);
            }

            if (mode == "gonzalez") {
                return GonzalezCenters(
                    G, k, delta, num_buckets, seed, single_core, sssp_verbose, sssp_timing);
            }

            std::cout << "ERROR: Unknown initializer '" << mode
                      << "'. Supported initializers are: random, greedy, gonzalez."
                      << std::endl;
            std::exit(-1);
        }

    }  // namespace kmedian
}  // namespace gbbs
