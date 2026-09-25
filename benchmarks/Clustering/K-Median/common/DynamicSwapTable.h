#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

#include "benchmarks/Clustering/K-Median/common/RangeMaxTree.h"
#include "benchmarks/Clustering/K-Median/common/StaticSwapTable.h"
#include "parlay/sequence.h"

namespace gbbs {
namespace kmedian {

class DynamicSwapTableTestAccess;

class DynamicSwapTable {
public:
    DynamicSwapTable() = default;

    template <class SwapTable>
    DynamicSwapTable(const SwapTable& static_table, size_t n) {
        Build(static_table, n);
    }

    template <class SwapTable>
    void Build(const SwapTable& static_table, size_t n) {
        num_rows_ = static_table.num_rows();
        num_cols_ = static_table.num_cols();

        if (num_rows_ == 0) {
            Fail("DynamicSwapTable requires at least one center.");
        }

        centers_ = parlay::sequence<uintE>(num_rows_);
        incoming_centers_ = parlay::sequence<uintE>(num_cols_);
        center_row_ = parlay::sequence<size_t>(n, kNoIndex);
        incoming_column_ = parlay::sequence<size_t>(n, kNoIndex);
        row_tie_break_keys_ = parlay::sequence<size_t>(num_rows_);

        for (size_t row = 0; row < num_rows_; ++row) {
            const uintE x = static_table.outgoing_center(row);
            if (x >= n) {
                Fail("Outgoing center is outside the vertex range.");
            }
            centers_[row] = x;
            center_row_[x] = row;
            row_tie_break_keys_[row] = static_cast<size_t>(x);
        }

        columns_.clear();
        columns_.reserve(num_cols_);
        column_maxima_.clear();

        for (size_t col = 0; col < num_cols_; ++col) {
            const uintE y = static_table.incoming_center(col);
            if (y >= n) {
                Fail("Incoming center is outside the vertex range.");
            }

            incoming_centers_[col] = y;
            incoming_column_[y] = col;

            parlay::sequence<double> values(num_rows_);
            for (size_t row = 0; row < num_rows_; ++row) {
                values[row] = static_table.gain(row, col);
            }

            columns_.emplace_back(
                values,
                row_tie_break_keys_.data());

            InsertColumnMaximum(col);
        }
    }

    size_t num_rows() const {
        return num_rows_;
    }

    size_t num_cols() const {
        return num_cols_;
    }

    uintE outgoing_center(size_t row) const {
        CheckRow(row);
        return centers_[row];
    }

    uintE incoming_center(size_t col) const {
        CheckColumn(col);
        return incoming_centers_[col];
    }

    size_t row_of_center(uintE x) const {
        return x < center_row_.size() ? center_row_[x] : kNoIndex;
    }

    size_t column_of_noncenter(uintE y) const {
        return y < incoming_column_.size() ? incoming_column_[y] : kNoIndex;
    }

    double gain(size_t row, size_t col) {
        CheckRow(row);
        CheckColumn(col);
        return columns_[col].Value(row);
    }

    void PointAdd(size_t col, size_t row, double delta) {
        CheckColumn(col);
        CheckRow(row);
        RemoveColumnMaximum(col);
        columns_[col].PointAdd(row, delta);
        InsertColumnMaximum(col);
    }

    void RangeAdd(size_t col, size_t left, size_t right, double delta) {
        CheckColumn(col);
        if (left > right) {
            return;
        }
        CheckRow(left);
        CheckRow(right);
        RemoveColumnMaximum(col);
        columns_[col].RangeAdd(left, right, delta);
        InsertColumnMaximum(col);
    }

    void RefreshColumnMaximum(size_t col) {
        CheckColumn(col);
        RemoveColumnMaximum(col);
        InsertColumnMaximum(col);
    }

    template <class OldDistance, class NewDistance>
    void ApplyDistanceDecrease(uintE p, uintE q, KMedianSolutionState& state,
                               OldDistance old_distance, NewDistance new_distance) {
        ValidateEndpoints(p, q, state);

        const bool p_is_center = state.is_center[p];
        const bool q_is_center = state.is_center[q];

        if (!p_is_center && !q_is_center) {
            ApplyNoncenterDistanceDecrease(p, q, state, old_distance, new_distance);
        } else if (p_is_center && q_is_center) {
            ApplyCenterDistanceDecrease(p, q, state, old_distance, new_distance);
        } else {
            ApplyMixedDistanceDecrease(p, q, state, old_distance, new_distance);
        }

        UpdateSolutionStateAfterDistanceDecrease(p, q, state, new_distance);
    }

    template <class OldDistance, class NewDistance>
    void ApplyNoncenterDistanceDecrease(uintE p, uintE q, const KMedianSolutionState& state,
                                        OldDistance old_distance, NewDistance new_distance) {
        ValidateEndpoints(p, q, state);
        if (state.is_center[p] || state.is_center[q]) {
            Fail("ApplyNoncenterDistanceDecrease requires two noncenters.");
        }

        const size_t affected_columns[2] = {column_of_noncenter(p), column_of_noncenter(q)};
        if (affected_columns[0] == kNoIndex || affected_columns[1] == kNoIndex) {
            Fail("Could not locate a noncenter column.");
        }

        for (const size_t col : affected_columns) {
            const uintE y = incoming_centers_[col];
            RemoveColumnMaximum(col);

            for (size_t row = 0; row < num_rows_; ++row) {
                const uintE x = centers_[row];
                const double diff_gain = GainContributionChange(
                    p, x, y, state, old_distance, new_distance)
                    + GainContributionChange(q, x, y, state, old_distance, new_distance);
                columns_[col].PointAdd(row, diff_gain);
            }

            InsertColumnMaximum(col);
        }
    }

    template <class OldDistance, class NewDistance>
    void ApplyCenterDistanceDecrease(uintE p, uintE q, const KMedianSolutionState& state,
                                     OldDistance old_distance, NewDistance new_distance) {
        ValidateEndpoints(p, q, state);
        if (!state.is_center[p] || !state.is_center[q]) {
            Fail("ApplyCenterDistanceDecrease requires two centers.");
        }

        const size_t p_row = row_of_center(p);
        const size_t q_row = row_of_center(q);
        if (p_row == kNoIndex || q_row == kNoIndex) {
            Fail("Could not locate a center row.");
        }

        for (size_t col = 0; col < num_cols_; ++col) {
            const uintE y = incoming_centers_[col];
            RemoveColumnMaximum(col);
            columns_[col].PointAdd(
                p_row, CenterRemovalContributionChange(p, q, y, state, old_distance, new_distance));
            columns_[col].PointAdd(
                q_row, CenterRemovalContributionChange(q, p, y, state, old_distance, new_distance));
            InsertColumnMaximum(col);
        }
    }

    template <class OldDistance, class NewDistance>
    void ApplyMixedDistanceDecrease(uintE p, uintE q, const KMedianSolutionState& state,
                                    OldDistance old_distance, NewDistance new_distance) {
        ValidateEndpoints(p, q, state);
        if (state.is_center[p] == state.is_center[q]) {
            Fail("ApplyMixedDistanceDecrease requires exactly one center endpoint.");
        }
        if (!state.is_center[p]) {
            std::swap(p, q);
        }

        const size_t p_row = row_of_center(p);
        const size_t q_col = column_of_noncenter(q);
        const uintE closest_center = state.nearest_center[q];
        const size_t closest_row = row_of_center(closest_center);
        if (p_row == kNoIndex || q_col == kNoIndex || closest_row == kNoIndex) {
            Fail("Could not locate mixed-update row, column, or closest center.");
        }

        const double old_q_to_closest = static_cast<double>(state.nearest_dist[q]);
        const double new_pq = static_cast<double>(new_distance(p, q));

        for (size_t col = 0; col < num_cols_; ++col) {
            const uintE y = incoming_centers_[col];
            RemoveColumnMaximum(col);

            if (col == q_col) {
                UpdateExplicitColumn(col, p, q, y, state, old_distance, new_distance);
                InsertColumnMaximum(col);
                continue;
            }

            const double p_row_diff = SinglePairGainContributionChange(
                p, q, p, y, state, old_distance, new_distance)
                + SinglePairGainContributionChange(q, p, p, y, state, old_distance, new_distance);
            columns_[col].PointAdd(p_row, p_row_diff);

            const double old_qy = static_cast<double>(old_distance(q, y));
            double common_diff_gain = 0.0;
            if (old_q_to_closest > new_pq) {
                common_diff_gain = new_pq <= old_qy
                    ? std::min(old_q_to_closest, old_qy) - old_q_to_closest
                    : new_pq - old_q_to_closest;
            }
            AddToRowsExcluding(col, p_row, closest_row, common_diff_gain);

            if (closest_center != p) {
                const double closest_diff = SinglePairGainContributionChange(
                    q, p, closest_center, y, state, old_distance, new_distance);
                columns_[col].PointAdd(closest_row, closest_diff);
            }

            InsertColumnMaximum(col);
        }
    }

    double max_gain() const {
        CheckBuilt();
        return column_maxima_.begin()->gain;
    }

    uintE best_outgoing_center() const {
        CheckBuilt();
        return column_maxima_.begin()->outgoing_center;
    }

    uintE best_incoming_center() const {
        CheckBuilt();
        return column_maxima_.begin()->incoming_center;
    }

private:
    friend class DynamicSwapTableTestAccess;

    struct ColumnMaximum {
        double gain = -std::numeric_limits<double>::infinity();
        uintE outgoing_center = kNoCenter;
        uintE incoming_center = kNoCenter;
        size_t column = kNoIndex;
    };

    struct BetterColumnMaximum {
        bool operator()(const ColumnMaximum& a, const ColumnMaximum& b) const {
            if (a.gain != b.gain) {
                return a.gain > b.gain;
            }
            if (a.outgoing_center != b.outgoing_center) {
                return a.outgoing_center < b.outgoing_center;
            }
            if (a.incoming_center != b.incoming_center) {
                return a.incoming_center < b.incoming_center;
            }
            return a.column < b.column;
        }
    };

    [[noreturn]] static void Fail(const char* message) {
        std::cout << "ERROR: " << message << std::endl;
        std::exit(-1);
    }

    void ValidateEndpoints(uintE p, uintE q, const KMedianSolutionState& state) const {
        if (p == q) {
            Fail("Distance decrease requires two distinct points.");
        }
        if (p >= state.is_center.size() || q >= state.is_center.size()) {
            Fail("Distance-decrease endpoint is outside the point set.");
        }
    }

    void CheckBuilt() const {
        if (num_rows_ == 0 || num_cols_ == 0 || column_maxima_.empty()) {
            Fail("DynamicSwapTable has not been initialized.");
        }
    }

    void CheckRow(size_t row) const {
        if (row >= num_rows_) {
            Fail("DynamicSwapTable row is outside the valid range.");
        }
    }

    void CheckColumn(size_t col) const {
        if (col >= num_cols_) {
            Fail("DynamicSwapTable column is outside the valid range.");
        }
    }

    template <class DistanceFunction>
    double DistanceToSwappedCenters(uintE a, uintE x, uintE y,
                                    const KMedianSolutionState& state,
                                    DistanceFunction distance) const {
        const double distance_to_y = static_cast<double>(distance(a, y));
        const Distance fallback = x == state.nearest_center[a]
            ? state.second_nearest_dist[a]
            : state.nearest_dist[a];
        return std::min(static_cast<double>(fallback), distance_to_y);
    }

    template <class OldDistance, class NewDistance>
    double GainContributionChange(uintE a, uintE x, uintE y,
                                  const KMedianSolutionState& state,
                                  OldDistance old_distance, NewDistance new_distance) const {
        const double distance_to_centers = static_cast<double>(state.nearest_dist[a]);
        const double old_after = DistanceToSwappedCenters(a, x, y, state, old_distance);
        const double new_after = DistanceToSwappedCenters(a, x, y, state, new_distance);
        return old_after - new_after;
    }

    template <class OldDistance, class NewDistance>
    double CenterRemovalContributionChange(uintE removed_center, uintE other_center, uintE incoming,
                                           const KMedianSolutionState& state,
                                           OldDistance old_distance, NewDistance new_distance) const {
        const double old_fallback = static_cast<double>(state.second_nearest_dist[removed_center]);
        const double new_fallback = std::min(
            old_fallback, static_cast<double>(new_distance(removed_center, other_center)));
        const double old_after = std::min(
            old_fallback, static_cast<double>(old_distance(removed_center, incoming)));
        const double new_after = std::min(
            new_fallback, static_cast<double>(new_distance(removed_center, incoming)));
        return old_after - new_after;
    }

    template <class OldDistance, class NewDistance>
    double SinglePairGainContributionChange(uintE a, uintE other, uintE x, uintE y,
                                            const KMedianSolutionState& state,
                                            OldDistance old_distance, NewDistance new_distance) const {
        const double old_to_centers = static_cast<double>(state.nearest_dist[a]);
        double new_to_centers = old_to_centers;
        if (state.is_center[other]) {
            new_to_centers = std::min(new_to_centers, static_cast<double>(new_distance(a, other)));
        }

        const double old_after = DistanceToSwappedCenters(a, x, y, state, old_distance);
        double new_after = old_after;
        const bool other_in_swapped_centers = (state.is_center[other] && other != x) || other == y;
        if (other_in_swapped_centers) {
            new_after = std::min(new_after, static_cast<double>(new_distance(a, other)));
        }

        return (new_to_centers - new_after) - (old_to_centers - old_after);
    }

    template <class OldDistance, class NewDistance>
    void UpdateExplicitColumn(size_t col, uintE p, uintE q, uintE y,
                              const KMedianSolutionState& state,
                              OldDistance old_distance, NewDistance new_distance) {
        for (size_t row = 0; row < num_rows_; ++row) {
            const uintE x = centers_[row];
            const double diff_gain = SinglePairGainContributionChange(
                p, q, x, y, state, old_distance, new_distance)
                + SinglePairGainContributionChange(q, p, x, y, state, old_distance, new_distance);
            columns_[col].PointAdd(row, diff_gain);
        }
    }

    void AddToRowsExcluding(size_t col, size_t first, size_t second, double delta) {
        if (delta == 0.0) {
            return;
        }
        if (first > second) {
            std::swap(first, second);
        }
        if (first > 0) {
            columns_[col].RangeAdd(0, first - 1, delta);
        }
        if (first + 1 < second) {
            columns_[col].RangeAdd(first + 1, second - 1, delta);
        }
        if (second + 1 < num_rows_) {
            columns_[col].RangeAdd(second + 1, num_rows_ - 1, delta);
        }
    }

    template <class NewDistance>
    void UpdateSolutionStateAfterDistanceDecrease(uintE p, uintE q, KMedianSolutionState& state,
                                                  NewDistance new_distance) const {
        const double old_p = static_cast<double>(state.nearest_dist[p]);
        const double old_q = static_cast<double>(state.nearest_dist[q]);
        UpdatePointSolutionState(p, q, state, new_distance);
        UpdatePointSolutionState(q, p, state, new_distance);
        state.cost += static_cast<double>(state.nearest_dist[p]) - old_p
            + static_cast<double>(state.nearest_dist[q]) - old_q;
    }

    template <class NewDistance>
    void UpdatePointSolutionState(uintE a, uintE other, KMedianSolutionState& state,
                                  NewDistance new_distance) const {
        if (!state.is_center[other]) {
            return;
        }

        const Distance candidate = static_cast<Distance>(new_distance(a, other));
        if (state.nearest_center[a] == other) {
            state.nearest_dist[a] = candidate;
            return;
        }
        if (better_center(candidate, other, state.nearest_dist[a], state.nearest_center[a])) {
            state.second_nearest_dist[a] = state.nearest_dist[a];
            state.second_nearest_center[a] = state.nearest_center[a];
            state.nearest_dist[a] = candidate;
            state.nearest_center[a] = other;
            return;
        }
        if (state.second_nearest_center[a] == other) {
            state.second_nearest_dist[a] = candidate;
            return;
        }
        if (better_center(candidate, other,
                          state.second_nearest_dist[a], state.second_nearest_center[a])) {
            state.second_nearest_dist[a] = candidate;
            state.second_nearest_center[a] = other;
        }
    }

    void RemoveColumnMaximum(size_t col) {
        const size_t row = columns_[col].MaxIndex();
        const ColumnMaximum maximum{
            columns_[col].MaxValue(), centers_[row], incoming_centers_[col], col};
        const auto it = column_maxima_.find(maximum);
        if (it == column_maxima_.end()) {
            Fail("Could not find old column maximum.");
        }
        column_maxima_.erase(it);
    }

    void InsertColumnMaximum(size_t col) {
        const size_t row = columns_[col].MaxIndex();
        column_maxima_.insert(ColumnMaximum{
            columns_[col].MaxValue(), centers_[row], incoming_centers_[col], col});
    }

    static constexpr size_t kNoIndex = std::numeric_limits<size_t>::max();

    size_t num_rows_ = 0;
    size_t num_cols_ = 0;
    parlay::sequence<uintE> centers_;
    parlay::sequence<uintE> incoming_centers_;
    parlay::sequence<size_t> center_row_;
    parlay::sequence<size_t> incoming_column_;
    parlay::sequence<size_t> row_tie_break_keys_;
    std::vector<RangeMaxTree> columns_;
    std::set<ColumnMaximum, BetterColumnMaximum> column_maxima_;
};

}  // namespace kmedian
}  // namespace gbbs
