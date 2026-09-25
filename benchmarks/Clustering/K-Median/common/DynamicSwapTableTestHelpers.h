#pragma once

#include <algorithm>
#include <limits>

#include "benchmarks/Clustering/K-Median/common/DynamicSwapTable.h"

namespace gbbs {
namespace kmedian {

class DynamicSwapTableTestAccess {
public:
    template <class OldDistance, class NewDistance>
    static void ApplyMixedDistanceDecreaseReference(
        DynamicSwapTable& table, uintE p, uintE q, const KMedianSolutionState& state,
        OldDistance old_distance, NewDistance new_distance) {

        if (p == q || p >= state.is_center.size() || q >= state.is_center.size()) {
            DynamicSwapTable::Fail("Invalid reference distance-decrease endpoints.");
        }
        if (state.is_center[p] == state.is_center[q]) {
            DynamicSwapTable::Fail("Reference mixed update requires exactly one center endpoint.");
        }
        if (!state.is_center[p]) {
            std::swap(p, q);
        }

        for (size_t col = 0; col < table.num_cols_; ++col) {
            const uintE y = table.incoming_centers_[col];
            table.RemoveColumnMaximum(col);

            for (size_t row = 0; row < table.num_rows_; ++row) {
                const uintE x = table.centers_[row];
                const double diff_gain = GainContributionChange(
                    table, p, x, y, old_distance, new_distance)
                    + GainContributionChange(table, q, x, y, old_distance, new_distance);
                table.columns_[col].PointAdd(row, diff_gain);
            }

            table.InsertColumnMaximum(col);
        }
    }

private:
    template <class DistanceFunction>
    static double DistanceToCenters(
        const DynamicSwapTable& table, uintE a, DistanceFunction distance) {

        double best = std::numeric_limits<double>::infinity();
        for (const uintE c : table.centers_) {
            best = std::min(best, static_cast<double>(distance(a, c)));
        }
        return best;
    }

    template <class DistanceFunction>
    static double DistanceToSwappedCenters(
        const DynamicSwapTable& table, uintE a, uintE x, uintE y, DistanceFunction distance) {

        double best = std::numeric_limits<double>::infinity();
        for (const uintE c : table.centers_) {
            if (c != x) {
                best = std::min(best, static_cast<double>(distance(a, c)));
            }
        }
        return std::min(best, static_cast<double>(distance(a, y)));
    }

    template <class OldDistance, class NewDistance>
    static double GainContributionChange(
        const DynamicSwapTable& table, uintE a, uintE x, uintE y,
        OldDistance old_distance, NewDistance new_distance) {

        const double old_to_centers = DistanceToCenters(table, a, old_distance);
        const double new_to_centers = DistanceToCenters(table, a, new_distance);
        const double old_after = DistanceToSwappedCenters(table, a, x, y, old_distance);
        const double new_after = DistanceToSwappedCenters(table, a, x, y, new_distance);
        return (new_to_centers - new_after) - (old_to_centers - old_after);
    }
};

}  // namespace kmedian
}  // namespace gbbs
