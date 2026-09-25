#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

#include "benchmarks/Clustering/K-Median/common/KMedianCommon.h"
#include "parlay/sequence.h"

namespace gbbs {
    namespace kmedian {

        struct SwapColumn {
            uintE incoming_center = kNoCenter;

            // gains[i] corresponds to removing centers[i].
            parlay::sequence<double> gains;

            double max_gain = -std::numeric_limits<double>::infinity();
            uintE best_outgoing_center = kNoCenter;
        };


        // Compute one swap-table column using an already available exact
        // distance metric. The metric must provide:
        //
        //     size_t n() const;
        //     distance(u, v)
        template <class Metric>
        SwapColumn ComputeSwapColumnFromDistances(
            const Metric& metric,
            const KMedianSolutionState& state,
            uintE y) {

            const size_t n = metric.n();
            const size_t k = state.centers.size();
            const Distance kInf = distance_infinity();

            if (y >= n) {
                std::cout << "ERROR: Incoming center " << y
                          << " is outside metric vertex range."
                          << std::endl;
                std::exit(-1);
            }

            if (state.is_center.size() != n) {
                std::cout << "ERROR: Solution state size does not match metric size."
                          << std::endl;
                std::exit(-1);
            }

            if (state.is_center[y]) {
                std::cout << "ERROR: Incoming center " << y
                          << " is already in C."
                          << std::endl;
                std::exit(-1);
            }

            double common_gain = 0.0;
            parlay::sequence<double> correction(k, 0.0);
            parlay::sequence<size_t> center_row(
                n,
                std::numeric_limits<size_t>::max());

            for (size_t i = 0; i < k; ++i) {
                center_row[state.centers[i]] = i;
            }

            for (size_t v = 0; v < n; ++v) {
                const Distance d1 = state.nearest_dist[v];
                const Distance d2 = state.second_nearest_dist[v];
                const auto raw_dy =
                    metric.distance(y, static_cast<uintE>(v));

                if (d1 == kInf || is_unreachable(raw_dy)) {
                    std::cout << "ERROR: Swap-table computation currently expects "
                                 "all vertices to be reachable."
                              << std::endl;
                    std::exit(-1);
                }

                const Distance dy =
                    static_cast<Distance>(raw_dy);

                const uintE nearest = state.nearest_center[v];
                const size_t row = center_row[nearest];

                if (row == std::numeric_limits<size_t>::max()) {
                    std::cout << "ERROR: Invalid nearest-center state."
                              << std::endl;
                    std::exit(-1);
                }

                const Distance common_new_dist =
                    std::min(d1, dy);

                const double common_contribution =
                    static_cast<double>(d1) -
                    static_cast<double>(common_new_dist);

                common_gain += common_contribution;

                const Distance special_new_dist =
                    std::min(d2, dy);

                const double special_contribution =
                    static_cast<double>(d1) -
                    static_cast<double>(special_new_dist);

                correction[row] +=
                    special_contribution - common_contribution;
            }

            SwapColumn column;
            column.incoming_center = y;
            column.gains = parlay::sequence<double>(k);

            for (size_t i = 0; i < k; ++i) {
                const double gain =
                    common_gain + correction[i];

                column.gains[i] = gain;

                const uintE x = state.centers[i];

                if (gain > column.max_gain ||
                    (gain == column.max_gain &&
                     x < column.best_outgoing_center)) {

                    column.max_gain = gain;
                    column.best_outgoing_center = x;
                }
            }

            return column;
        }

        class StaticSwapTable {
            public:
                StaticSwapTable() = default;

                template <class Metric>
                StaticSwapTable(
                    const Metric& metric,
                    const KMedianSolutionState& state) {

                    BuildFromDistances(metric, state);
                }


                template <class Metric>
                void BuildFromDistances(
                    const Metric& metric,
                    const KMedianSolutionState& state) {

                    const size_t n = metric.n();

                    if (state.is_center.size() != n) {
                        std::cout << "ERROR: Solution state size does not match metric size."
                                  << std::endl;
                        std::exit(-1);
                    }

                    centers_ = state.centers;
                    incoming_centers_.clear();

                    for (uintE v = 0; v < n; ++v) {
                        if (!state.is_center[v]) {
                            incoming_centers_.push_back(v);
                        }
                    }

                    num_rows_ = centers_.size();
                    num_cols_ = incoming_centers_.size();

                    gains_ = parlay::sequence<double>(
                        num_rows_ * num_cols_,
                        -std::numeric_limits<double>::infinity());

                    max_gain_ =
                        -std::numeric_limits<double>::infinity();
                    best_outgoing_center_ = kNoCenter;
                    best_incoming_center_ = kNoCenter;

                    for (size_t col = 0; col < num_cols_; ++col) {
                        const uintE y =
                            incoming_centers_[col];

                        auto column =
                            ComputeSwapColumnFromDistances(
                                metric,
                                state,
                                y);

                        for (size_t row = 0; row < num_rows_; ++row) {
                            gains_[index(row, col)] =
                                column.gains[row];
                        }

                        const double gain =
                            column.max_gain;
                        const uintE x =
                            column.best_outgoing_center;

                        if (gain > max_gain_ ||
                            (gain == max_gain_ &&
                             better_swap(
                                 x,
                                 y,
                                 best_outgoing_center_,
                                 best_incoming_center_))) {

                            max_gain_ = gain;
                            best_outgoing_center_ = x;
                            best_incoming_center_ = y;
                        }
                    }
                }

                size_t num_rows() const {
                    return num_rows_;
                }

                size_t num_cols() const {
                    return num_cols_;
                }

                double gain(size_t row, size_t col) const {
                    return gains_[index(row, col)];
                }

                uintE outgoing_center(size_t row) const {
                    return centers_[row];
                }

                uintE incoming_center(size_t col) const {
                    return incoming_centers_[col];
                }

                double max_gain() const {
                    return max_gain_;
                }

                uintE best_outgoing_center() const {
                    return best_outgoing_center_;
                }

                uintE best_incoming_center() const {
                    return best_incoming_center_;
                }

            private:
                size_t index(size_t row, size_t col) const {
                    return row * num_cols_ + col;
                }

                static bool better_swap(
                    uintE candidate_x,
                    uintE candidate_y,
                    uintE current_x,
                    uintE current_y) {

                    if (current_x == kNoCenter) {
                        return true;
                    }
                    if (candidate_x < current_x) {
                        return true;
                    }
                    if (candidate_x > current_x) {
                        return false;
                    }

                    return candidate_y < current_y;
                }

                size_t num_rows_ = 0;
                size_t num_cols_ = 0;

                parlay::sequence<uintE> centers_;
                parlay::sequence<uintE> incoming_centers_;

                // Row-major:
                //
                // gains_[row * num_cols_ + col]
                parlay::sequence<double> gains_;

                double max_gain_ = -std::numeric_limits<double>::infinity();
                uintE best_outgoing_center_ = kNoCenter;
                uintE best_incoming_center_ = kNoCenter;
        };

        template <class Metric>
        void VerifySwapColumnBruteForce(
            const Metric& metric,
            const KMedianSolutionState& state,
            const SwapColumn& column) {

            const uintE y = column.incoming_center;
            const size_t k = state.centers.size();

            std::cout << "### Swap-column verification for y=" << y << std::endl;

            for (size_t row = 0; row < k; ++row) {
                const uintE x = state.centers[row];
                parlay::sequence<uintE> swapped_centers = state.centers;
                swapped_centers[row] = y;

                auto swapped_state = BuildSolutionStateFromDistances(
                    metric,
                    swapped_centers);

                const double brute_gain = state.cost - swapped_state.cost;
                const double fast_gain = column.gains[row];
                const double error = std::abs(fast_gain - brute_gain);

                std::cout << "remove=" << x
                          << " add=" << y
                          << " fast_gain=" << fast_gain
                          << " brute_gain=" << brute_gain
                          << " error=" << error
                          << std::endl;
            }
        }

        template <class Metric>
        void VerifyStaticSwapTableBruteForce(
            const Metric& metric,
            const KMedianSolutionState& state,
            const StaticSwapTable& table) {

            std::cout << "### Full swap-table brute-force verification" << std::endl;

            for (size_t row = 0; row < table.num_rows(); ++row) {
                const uintE x = table.outgoing_center(row);

                for (size_t col = 0; col < table.num_cols(); ++col) {
                    const uintE y = table.incoming_center(col);
                    parlay::sequence<uintE> swapped_centers = state.centers;
                    swapped_centers[row] = y;

                    auto swapped_state = BuildSolutionStateFromDistances(
                        metric,
                        swapped_centers);

                    const double brute_gain = state.cost - swapped_state.cost;
                    const double fast_gain = table.gain(row, col);
                    const double error = std::abs(fast_gain - brute_gain);

                    std::cout << "remove=" << x
                              << " add=" << y
                              << " fast_gain=" << fast_gain
                              << " brute_gain=" << brute_gain
                              << " error=" << error
                              << std::endl;
                }
            }
        }

    }  // namespace kmedian
}  // namespace gbbs
