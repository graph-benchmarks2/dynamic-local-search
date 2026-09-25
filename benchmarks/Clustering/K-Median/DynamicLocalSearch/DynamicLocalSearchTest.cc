#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include "benchmarks/Clustering/K-Median/DynamicLocalSearch/DynamicLocalSearch.h"
#include "benchmarks/Clustering/K-Median/common/KMedianCommon.h"
#include "parlay/sequence.h"

namespace {

using DistanceMatrix = std::vector<std::vector<double>>;

void CheckEqual(double actual, double expected, const char* label) {
    constexpr double kTolerance = 1e-7;

    if (std::abs(actual - expected) > kTolerance) {
        std::cout << "ERROR: " << label
                  << ": expected " << expected
                  << ", got " << actual
                  << std::endl;
        std::exit(-1);
    }
}

void CheckEqual(size_t actual, size_t expected, const char* label) {
    if (actual != expected) {
        std::cout << "ERROR: " << label
                  << ": expected " << expected
                  << ", got " << actual
                  << std::endl;
        std::exit(-1);
    }
}

void CheckEqual(gbbs::uintE actual, gbbs::uintE expected, const char* label) {
    if (actual != expected) {
        std::cout << "ERROR: " << label
                  << ": expected " << expected
                  << ", got " << actual
                  << std::endl;
        std::exit(-1);
    }
}

DistanceMatrix MakePathDistances(size_t n) {
    DistanceMatrix distances(n, std::vector<double>(n, 0.0));

    for (size_t u = 0; u < n; ++u) {
        for (size_t v = 0; v < n; ++v) {
            distances[u][v] = static_cast<double>(
                u > v ? u - v : v - u);
        }
    }

    return distances;
}

gbbs::kmedian::KMedianSolutionState BuildMetricSolutionState(
    const DistanceMatrix& distances,
    const parlay::sequence<gbbs::uintE>& centers) {

    const size_t n = distances.size();
    const auto kInf = gbbs::kmedian::distance_infinity();

    gbbs::kmedian::KMedianSolutionState state;
    state.centers = centers;
    state.is_center = parlay::sequence<bool>(n, false);
    state.nearest_center =
        parlay::sequence<gbbs::uintE>(n, gbbs::kmedian::kNoCenter);
    state.second_nearest_center =
        parlay::sequence<gbbs::uintE>(n, gbbs::kmedian::kNoCenter);
    state.nearest_dist =
        parlay::sequence<gbbs::kmedian::Distance>(n, kInf);
    state.second_nearest_dist =
        parlay::sequence<gbbs::kmedian::Distance>(n, kInf);

    for (const auto c : centers) {
        state.is_center[c] = true;
    }

    for (size_t v = 0; v < n; ++v) {
        for (const auto c : centers) {
            const auto d = static_cast<gbbs::kmedian::Distance>(
                distances[v][c]);

            if (gbbs::kmedian::better_center(
                    d,
                    c,
                    state.nearest_dist[v],
                    state.nearest_center[v])) {

                state.second_nearest_dist[v] =
                    state.nearest_dist[v];

                state.second_nearest_center[v] =
                    state.nearest_center[v];

                state.nearest_dist[v] = d;
                state.nearest_center[v] = c;
            } else if (gbbs::kmedian::better_center(
                           d,
                           c,
                           state.second_nearest_dist[v],
                           state.second_nearest_center[v])) {

                state.second_nearest_dist[v] = d;
                state.second_nearest_center[v] = c;
            }
        }

        state.cost += static_cast<double>(state.nearest_dist[v]);
    }

    return state;
}

class MetricSwapTable {
public:
    MetricSwapTable() = default;

    MetricSwapTable(
        const DistanceMatrix& distances,
        const gbbs::kmedian::KMedianSolutionState& state) {

        const size_t n = distances.size();

        centers_ = state.centers;

        for (gbbs::uintE y = 0; y < n; ++y) {
            if (!state.is_center[y]) {
                incoming_centers_.push_back(y);
            }
        }

        gains_ = parlay::sequence<double>(
            centers_.size() * incoming_centers_.size(),
            0.0);

        max_gain_ = -std::numeric_limits<double>::infinity();

        for (size_t row = 0; row < centers_.size(); ++row) {
            const auto x = centers_[row];

            for (size_t col = 0; col < incoming_centers_.size(); ++col) {
                const auto y = incoming_centers_[col];

                double swapped_cost = 0.0;

                for (size_t v = 0; v < n; ++v) {
                    double best = std::numeric_limits<double>::infinity();

                    for (const auto c : centers_) {
                        if (c != x) {
                            best = std::min(best, distances[v][c]);
                        }
                    }

                    best = std::min(best, distances[v][y]);
                    swapped_cost += best;
                }

                const double gain = state.cost - swapped_cost;
                gains_[Index(row, col)] = gain;

                if (gain > max_gain_
                    || (gain == max_gain_
                        && BetterSwap(
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
    }

    size_t num_rows() const {
        return centers_.size();
    }

    size_t num_cols() const {
        return incoming_centers_.size();
    }

    double gain(size_t row, size_t col) const {
        return gains_[Index(row, col)];
    }

    gbbs::uintE outgoing_center(size_t row) const {
        return centers_[row];
    }

    gbbs::uintE incoming_center(size_t col) const {
        return incoming_centers_[col];
    }

    double max_gain() const {
        return max_gain_;
    }

    gbbs::uintE best_outgoing_center() const {
        return best_outgoing_center_;
    }

    gbbs::uintE best_incoming_center() const {
        return best_incoming_center_;
    }

private:
    size_t Index(size_t row, size_t col) const {
        return row * incoming_centers_.size() + col;
    }

    static bool BetterSwap(
        gbbs::uintE candidate_x,
        gbbs::uintE candidate_y,
        gbbs::uintE current_x,
        gbbs::uintE current_y) {

        if (current_x == gbbs::kmedian::kNoCenter) {
            return true;
        }

        if (candidate_x != current_x) {
            return candidate_x < current_x;
        }

        return candidate_y < current_y;
    }

    parlay::sequence<gbbs::uintE> centers_;
    parlay::sequence<gbbs::uintE> incoming_centers_;
    parlay::sequence<double> gains_;

    double max_gain_ = -std::numeric_limits<double>::infinity();
    gbbs::uintE best_outgoing_center_ = gbbs::kmedian::kNoCenter;
    gbbs::uintE best_incoming_center_ = gbbs::kmedian::kNoCenter;
};

struct MetricRecoveryResult {
    gbbs::kmedian::KMedianSolutionState state;
    MetricSwapTable final_table;

    size_t num_swaps = 0;

    double initial_cost = std::numeric_limits<double>::infinity();
    double final_cost = std::numeric_limits<double>::infinity();
};

MetricRecoveryResult RunMetricStaticRecovery(
    const DistanceMatrix& distances,
    const parlay::sequence<gbbs::uintE>& initial_centers,
    double ls_delta) {

    MetricRecoveryResult result;

    auto state = BuildMetricSolutionState(
        distances,
        initial_centers);

    result.initial_cost = state.cost;

    while (true) {
        MetricSwapTable table(distances, state);

        const double threshold =
            (ls_delta / static_cast<double>(state.centers.size()))
            * state.cost;

        if (!(table.max_gain() > 0.0
              && table.max_gain() >= threshold)) {

            break;
        }

        auto centers = gbbs::kmedian::CentersAfterSwap(
            state.centers,
            table.best_outgoing_center(),
            table.best_incoming_center());

        state = BuildMetricSolutionState(
            distances,
            centers);

        ++result.num_swaps;
    }

    result.final_cost = state.cost;
    result.final_table = MetricSwapTable(distances, state);
    result.state = std::move(state);

    return result;
}

void TestNoSwapBatch() {
    constexpr double kDelta = 0.1;

    DistanceMatrix distances = MakePathDistances(8);
    parlay::sequence<gbbs::uintE> centers = {1, 5};

    auto state = BuildMetricSolutionState(
        distances,
        centers);

    MetricSwapTable initial_table(distances, state);
    gbbs::kmedian::DynamicSwapTable sw_decr(
        initial_table,
        distances.size());

    parlay::sequence<gbbs::kmedian::DistanceDecrease> batch = {
        gbbs::kmedian::DistanceDecrease{0, 7, 6.5}
    };

    auto get_distance = [&](gbbs::uintE u, gbbs::uintE v) {
        return distances[u][v];
    };

    auto set_distance =
        [&](gbbs::uintE u, gbbs::uintE v, double value) {

            distances[u][v] = value;
            distances[v][u] = value;
        };

    auto static_recovery =
        [&](const parlay::sequence<gbbs::uintE>& new_centers) {

            return RunMetricStaticRecovery(
                distances,
                new_centers,
                kDelta);
        };

    auto result = gbbs::kmedian::ProcessDistanceDecreaseBatch(
        sw_decr,
        state,
        batch,
        kDelta,
        get_distance,
        set_distance,
        static_recovery);

    CheckEqual(result.cost_before_batch, 8.0, "no-swap initial cost");
    CheckEqual(result.cost_after_batch_before_ls, 8.0, "no-swap post-batch cost");
    CheckEqual(result.final_cost, 8.0, "no-swap final cost");
    CheckEqual(result.num_swaps, size_t{0}, "no-swap number of swaps");

    if (result.triggered_swap) {
        std::cout << "ERROR: No-swap batch unexpectedly triggered a swap."
                  << std::endl;
        std::exit(-1);
    }

    std::cout << "DynamicLocalSearch no-swap batch test passed."
              << std::endl;
}

void TestSwapTriggeringBatchAndRebuild() {
    constexpr double kDelta = 0.1;

    DistanceMatrix distances = MakePathDistances(8);
    parlay::sequence<gbbs::uintE> centers = {1, 5};

    auto state = BuildMetricSolutionState(
        distances,
        centers);

    MetricSwapTable initial_table(distances, state);
    gbbs::kmedian::DynamicSwapTable sw_decr(
        initial_table,
        distances.size());

    parlay::sequence<gbbs::kmedian::DistanceDecrease> first_batch = {
        gbbs::kmedian::DistanceDecrease{0, 2, 1.0}
    };

    auto get_distance = [&](gbbs::uintE u, gbbs::uintE v) {
        return distances[u][v];
    };

    auto set_distance =
        [&](gbbs::uintE u, gbbs::uintE v, double value) {

            distances[u][v] = value;
            distances[v][u] = value;
        };

    auto static_recovery =
        [&](const parlay::sequence<gbbs::uintE>& new_centers) {

            return RunMetricStaticRecovery(
                distances,
                new_centers,
                kDelta);
        };

    auto first_result = gbbs::kmedian::ProcessDistanceDecreaseBatch(
        sw_decr,
        state,
        first_batch,
        kDelta,
        get_distance,
        set_distance,
        static_recovery);

    if (!first_result.triggered_swap) {
        std::cout << "ERROR: Swap-triggering batch did not trigger."
                  << std::endl;
        std::exit(-1);
    }

    CheckEqual(
        first_result.first_outgoing_center,
        gbbs::uintE{1},
        "first dynamic outgoing center");

    CheckEqual(
        first_result.first_incoming_center,
        gbbs::uintE{2},
        "first dynamic incoming center");

    CheckEqual(
        first_result.trigger_max_gain,
        1.0,
        "first dynamic gain");

    CheckEqual(
        first_result.cost_after_first_swap,
        7.0,
        "cost after first dynamic swap");

    CheckEqual(
        first_result.final_cost,
        7.0,
        "first batch final cost");

    CheckEqual(
        first_result.num_swaps,
        size_t{1},
        "first batch total swaps");

    CheckEqual(
        first_result.num_static_swaps,
        size_t{0},
        "first batch static recovery swaps");

    CheckEqual(
        state.centers[0],
        gbbs::uintE{2},
        "rebuilt state first center");

    CheckEqual(
        state.centers[1],
        gbbs::uintE{5},
        "rebuilt state second center");

    MetricSwapTable expected_after_rebuild(
        distances,
        state);

    CheckEqual(
        sw_decr.max_gain(),
        expected_after_rebuild.max_gain(),
        "rebuilt SW_decr max gain");

    // Process another batch using the rebuilt SW_decr: the second batch starts from C = {2, 5}.
    parlay::sequence<gbbs::kmedian::DistanceDecrease> second_batch = {
        gbbs::kmedian::DistanceDecrease{1, 6, 4.5}
    };

    auto second_result = gbbs::kmedian::ProcessDistanceDecreaseBatch(
        sw_decr,
        state,
        second_batch,
        kDelta,
        get_distance,
        set_distance,
        static_recovery);

    MetricSwapTable expected_second_table(
        distances,
        state);

    CheckEqual(
        sw_decr.max_gain(),
        expected_second_table.max_gain(),
        "second batch SW_decr max gain");

    CheckEqual(
        second_result.final_cost,
        state.cost,
        "second batch final cost");

    std::cout << "DynamicLocalSearch swap/rebuild/carry-over test passed."
              << std::endl;
}

void TestDynamicFirstSwapThenStaticRecoverySwap() {
    constexpr double kDelta = 0.1;

    DistanceMatrix distances = MakePathDistances(8);
    parlay::sequence<gbbs::uintE> centers = {1, 5};

    auto state = BuildMetricSolutionState(
        distances,
        centers);

    MetricSwapTable initial_table(distances, state);

    gbbs::kmedian::DynamicSwapTable sw_decr(
        initial_table,
        distances.size());

    parlay::sequence<gbbs::kmedian::DistanceDecrease> batch = {
        gbbs::kmedian::DistanceDecrease{2, 4, 0.5}
    };

    auto get_distance = [&](gbbs::uintE u, gbbs::uintE v) {
        return distances[u][v];
    };

    auto set_distance =
        [&](gbbs::uintE u, gbbs::uintE v, double value) {

            distances[u][v] = value;
            distances[v][u] = value;
        };

    auto static_recovery =
        [&](const parlay::sequence<gbbs::uintE>& new_centers) {

            return RunMetricStaticRecovery(
                distances,
                new_centers,
                kDelta);
        };

    auto result = gbbs::kmedian::ProcessDistanceDecreaseBatch(
        sw_decr,
        state,
        batch,
        kDelta,
        get_distance,
        set_distance,
        static_recovery);

    if (!result.triggered_swap) {
        std::cout << "ERROR: Multi-swap batch did not trigger a dynamic swap."
                  << std::endl;
        std::exit(-1);
    }

    CheckEqual(
        result.cost_before_batch,
        8.0,
        "multi-swap cost before batch");

    CheckEqual(
        result.cost_after_batch_before_ls,
        8.0,
        "multi-swap cost after batch before LS");

    CheckEqual(
        result.trigger_max_gain,
        0.5,
        "multi-swap dynamic first gain");

    CheckEqual(
        result.first_outgoing_center,
        gbbs::uintE{1},
        "multi-swap first outgoing center");

    CheckEqual(
        result.first_incoming_center,
        gbbs::uintE{2},
        "multi-swap first incoming center");

    CheckEqual(
        result.cost_after_first_swap,
        7.5,
        "multi-swap cost after dynamic first swap");

    CheckEqual(
        result.num_static_swaps,
        size_t{1},
        "multi-swap number of static recovery swaps");

    CheckEqual(
        result.num_swaps,
        size_t{2},
        "multi-swap total number of swaps");

    CheckEqual(
        result.final_cost,
        6.5,
        "multi-swap final cost");

    CheckEqual(
        state.centers[0],
        gbbs::uintE{2},
        "multi-swap final first center");

    CheckEqual(
        state.centers[1],
        gbbs::uintE{6},
        "multi-swap final second center");

    MetricSwapTable expected_final_table(
        distances,
        state);

    CheckEqual(
        sw_decr.max_gain(),
        expected_final_table.max_gain(),
        "multi-swap rebuilt SW_decr max gain");

    CheckEqual(
        result.final_max_gain,
        expected_final_table.max_gain(),
        "multi-swap final reported max gain");

    std::cout
        << "DynamicLocalSearch dynamic-first/static-recovery test passed."
        << std::endl;
}

}  // namespace

int main() {
    TestNoSwapBatch();
    TestSwapTriggeringBatchAndRebuild();
    TestDynamicFirstSwapThenStaticRecoverySwap();

    return 0;
}
