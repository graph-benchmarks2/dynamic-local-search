#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "gbbs/gbbs.h"

#include "benchmarks/Clustering/K-Median/common/DynamicSwapTable.h"
#include "benchmarks/Clustering/K-Median/common/DynamicSwapTableTestHelpers.h"
#include "benchmarks/Clustering/K-Median/common/KMedianCommon.h"
#include "benchmarks/Clustering/K-Median/common/StaticSwapTable.h"
#include "parlay/sequence.h"

namespace {

using DistanceMatrix = std::vector<std::vector<double>>;

void CheckEqual(double actual, double expected, const char* label) {
    constexpr double kTolerance = 1e-7;
    if (std::abs(actual - expected) > kTolerance) {
        std::cout << std::setprecision(17) << "ERROR: " << label
                  << ": expected " << expected << ", got " << actual
                  << ", abs_error = " << std::abs(actual - expected) << std::endl;
        std::exit(-1);
    }
}

void CheckEqual(gbbs::uintE actual, gbbs::uintE expected, const char* label) {
    if (actual != expected) {
        std::cout << "ERROR: " << label << ": expected " << expected
                  << ", got " << actual << std::endl;
        std::exit(-1);
    }
}

void CheckEqual(size_t actual, size_t expected, const char* label) {
    if (actual != expected) {
        std::cout << "ERROR: " << label << ": expected " << expected
                  << ", got " << actual << std::endl;
        std::exit(-1);
    }
}

DistanceMatrix MakePathDistances(size_t n) {
    DistanceMatrix distances(n, std::vector<double>(n, 0.0));
    for (size_t u = 0; u < n; ++u) {
        for (size_t v = 0; v < n; ++v) {
            distances[u][v] = static_cast<double>(u > v ? u - v : v - u);
        }
    }
    return distances;
}

struct MatrixMetric {
    const DistanceMatrix& distances;

    size_t n() const {
        return distances.size();
    }

    double distance(gbbs::uintE u, gbbs::uintE v) const {
        return distances[u][v];
    }
};

struct BruteBestSwap {
    double gain;
    size_t row;
    size_t col;
};

BruteBestSwap FindBruteBestSwap(
    const std::vector<std::vector<double>>& gains,
    const gbbs::kmedian::DynamicSwapTable& table) {

    BruteBestSwap best{-std::numeric_limits<double>::infinity(), 0, 0};
    bool initialized = false;

    for (size_t row = 0; row < table.num_rows(); ++row) {
        for (size_t col = 0; col < table.num_cols(); ++col) {
            const double gain = gains[row][col];
            if (!initialized) {
                best = {gain, row, col};
                initialized = true;
                continue;
            }

            const auto candidate = std::make_pair(
                table.outgoing_center(row), table.incoming_center(col));
            const auto current = std::make_pair(
                table.outgoing_center(best.row), table.incoming_center(best.col));
            if (gain > best.gain || (gain == best.gain && candidate < current)) {
                best = {gain, row, col};
            }
        }
    }
    return best;
}

void VerifyWholeTable(
    gbbs::kmedian::DynamicSwapTable& table,
    const std::vector<std::vector<double>>& brute_gains) {

    for (size_t row = 0; row < table.num_rows(); ++row) {
        for (size_t col = 0; col < table.num_cols(); ++col) {
            CheckEqual(table.gain(row, col), brute_gains[row][col], "dynamic swap-table entry");
        }
    }

    const BruteBestSwap best = FindBruteBestSwap(brute_gains, table);
    CheckEqual(table.max_gain(), best.gain, "dynamic global maximum gain");
    CheckEqual(table.best_outgoing_center(), table.outgoing_center(best.row), "best outgoing center");
    CheckEqual(table.best_incoming_center(), table.incoming_center(best.col), "best incoming center");
}

std::vector<std::vector<double>> ComputeBruteSwapTable(
    const DistanceMatrix& distances,
    const parlay::sequence<gbbs::uintE>& centers,
    const gbbs::kmedian::DynamicSwapTable& table) {

    std::vector<std::vector<double>> gains(
        table.num_rows(), std::vector<double>(table.num_cols(), 0.0));

    double current_cost = 0.0;
    for (size_t v = 0; v < distances.size(); ++v) {
        double best = std::numeric_limits<double>::infinity();
        for (const gbbs::uintE c : centers) {
            best = std::min(best, distances[v][c]);
        }
        current_cost += best;
    }

    for (size_t row = 0; row < table.num_rows(); ++row) {
        const gbbs::uintE x = table.outgoing_center(row);
        for (size_t col = 0; col < table.num_cols(); ++col) {
            const gbbs::uintE y = table.incoming_center(col);
            double swapped_cost = 0.0;

            for (size_t v = 0; v < distances.size(); ++v) {
                double best = std::numeric_limits<double>::infinity();
                for (const gbbs::uintE c : centers) {
                    if (c != x) {
                        best = std::min(best, distances[v][c]);
                    }
                }
                swapped_cost += std::min(best, distances[v][y]);
            }
            gains[row][col] = current_cost - swapped_cost;
        }
    }
    return gains;
}

void VerifyInitialization(
    gbbs::kmedian::DynamicSwapTable& dynamic_table,
    const gbbs::kmedian::StaticSwapTable& static_table) {

    CheckEqual(dynamic_table.num_rows(), static_table.num_rows(), "number of rows");
    CheckEqual(dynamic_table.num_cols(), static_table.num_cols(), "number of columns");

    std::vector<std::vector<double>> gains(
        static_table.num_rows(), std::vector<double>(static_table.num_cols()));

    for (size_t row = 0; row < static_table.num_rows(); ++row) {
        const auto x = static_table.outgoing_center(row);
        CheckEqual(dynamic_table.outgoing_center(row), x, "outgoing center");
        CheckEqual(dynamic_table.row_of_center(x), row, "center-to-row mapping");
        for (size_t col = 0; col < static_table.num_cols(); ++col) {
            gains[row][col] = static_table.gain(row, col);
        }
    }

    for (size_t col = 0; col < static_table.num_cols(); ++col) {
        const auto y = static_table.incoming_center(col);
        CheckEqual(dynamic_table.incoming_center(col), y, "incoming center");
        CheckEqual(dynamic_table.column_of_noncenter(y), col, "noncenter-to-column mapping");
    }

    VerifyWholeTable(dynamic_table, gains);
}

void RunBasicDistanceDecreaseTests() {
    const DistanceMatrix old_distances = MakePathDistances(6);
    const MatrixMetric metric{old_distances};
    parlay::sequence<gbbs::uintE> centers = {1, 4};
    auto state = gbbs::kmedian::BuildSolutionStateFromDistances(metric, centers);
    gbbs::kmedian::StaticSwapTable static_table(metric, state);
    auto old_distance = [&](gbbs::uintE u, gbbs::uintE v) { return old_distances[u][v]; };

    gbbs::kmedian::DynamicSwapTable init_table(static_table, old_distances.size());
    VerifyInitialization(init_table, static_table);
    std::cout << "DynamicSwapTable initialization test passed." << std::endl;

    DistanceMatrix noncenter_distances = old_distances;
    noncenter_distances[0][5] = noncenter_distances[5][0] = 1.0;
    auto noncenter_new = [&](gbbs::uintE u, gbbs::uintE v) { return noncenter_distances[u][v]; };
    gbbs::kmedian::DynamicSwapTable noncenter_table(static_table, old_distances.size());
    noncenter_table.ApplyNoncenterDistanceDecrease(0, 5, state, old_distance, noncenter_new);
    VerifyWholeTable(noncenter_table, ComputeBruteSwapTable(noncenter_distances, centers, noncenter_table));
    std::cout << "DynamicSwapTable noncenter/noncenter decrease test passed." << std::endl;

    DistanceMatrix center_distances = old_distances;
    center_distances[1][4] = center_distances[4][1] = 1.0;
    auto center_new = [&](gbbs::uintE u, gbbs::uintE v) { return center_distances[u][v]; };
    gbbs::kmedian::DynamicSwapTable center_table(static_table, old_distances.size());
    center_table.ApplyCenterDistanceDecrease(1, 4, state, old_distance, center_new);
    VerifyWholeTable(center_table, ComputeBruteSwapTable(center_distances, centers, center_table));
    std::cout << "DynamicSwapTable center/center decrease test passed." << std::endl;

    DistanceMatrix mixed_distances = old_distances;
    mixed_distances[1][5] = mixed_distances[5][1] = 1.0;
    auto mixed_new = [&](gbbs::uintE u, gbbs::uintE v) { return mixed_distances[u][v]; };
    const auto mixed_expected = ComputeBruteSwapTable(mixed_distances, centers, init_table);

    gbbs::kmedian::DynamicSwapTable mixed_reference(static_table, old_distances.size());
    gbbs::kmedian::DynamicSwapTableTestAccess::ApplyMixedDistanceDecreaseReference(
        mixed_reference, 1, 5, state, old_distance, mixed_new);
    VerifyWholeTable(mixed_reference, mixed_expected);
    std::cout << "DynamicSwapTable mixed reference decrease test passed." << std::endl;

    gbbs::kmedian::DynamicSwapTable reversed_reference(static_table, old_distances.size());
    gbbs::kmedian::DynamicSwapTableTestAccess::ApplyMixedDistanceDecreaseReference(
        reversed_reference, 5, 1, state, old_distance, mixed_new);
    VerifyWholeTable(reversed_reference, mixed_expected);
    std::cout << "DynamicSwapTable mixed reversed-endpoint test passed." << std::endl;
}

void RunOptimizedMixedDistanceDecreaseTests() {
    const DistanceMatrix old_distances = MakePathDistances(8);
    const MatrixMetric metric{old_distances};
    parlay::sequence<gbbs::uintE> centers = {1, 4, 6};
    auto state = gbbs::kmedian::BuildSolutionStateFromDistances(metric, centers);
    gbbs::kmedian::StaticSwapTable static_table(metric, state);
    auto old_distance = [&](gbbs::uintE u, gbbs::uintE v) { return old_distances[u][v]; };

    DistanceMatrix mixed = old_distances;
    mixed[1][7] = mixed[7][1] = 1.0;
    auto mixed_new = [&](gbbs::uintE u, gbbs::uintE v) { return mixed[u][v]; };
    gbbs::kmedian::DynamicSwapTable optimized(static_table, old_distances.size());
    gbbs::kmedian::DynamicSwapTable reference(static_table, old_distances.size());
    optimized.ApplyMixedDistanceDecrease(1, 7, state, old_distance, mixed_new);
    gbbs::kmedian::DynamicSwapTableTestAccess::ApplyMixedDistanceDecreaseReference(
        reference, 1, 7, state, old_distance, mixed_new);
    const auto expected = ComputeBruteSwapTable(mixed, centers, optimized);
    VerifyWholeTable(reference, expected);
    VerifyWholeTable(optimized, expected);
    std::cout << "DynamicSwapTable optimized mixed decrease test passed." << std::endl;

    gbbs::kmedian::DynamicSwapTable reversed(static_table, old_distances.size());
    reversed.ApplyMixedDistanceDecrease(7, 1, state, old_distance, mixed_new);
    VerifyWholeTable(reversed, expected);
    std::cout << "DynamicSwapTable optimized mixed reversed-endpoint test passed." << std::endl;

    DistanceMatrix closest = old_distances;
    closest[6][7] = closest[7][6] = 0.25;
    auto closest_new = [&](gbbs::uintE u, gbbs::uintE v) { return closest[u][v]; };
    gbbs::kmedian::DynamicSwapTable closest_opt(static_table, old_distances.size());
    gbbs::kmedian::DynamicSwapTable closest_ref(static_table, old_distances.size());
    closest_opt.ApplyMixedDistanceDecrease(6, 7, state, old_distance, closest_new);
    gbbs::kmedian::DynamicSwapTableTestAccess::ApplyMixedDistanceDecreaseReference(
        closest_ref, 6, 7, state, old_distance, closest_new);
    const auto closest_expected = ComputeBruteSwapTable(closest, centers, closest_opt);
    VerifyWholeTable(closest_ref, closest_expected);
    VerifyWholeTable(closest_opt, closest_expected);
    std::cout << "DynamicSwapTable optimized mixed closest-center test passed." << std::endl;
}

void RunSequentialDistanceDecreaseTest() {
    DistanceMatrix distances = MakePathDistances(8);
    const MatrixMetric metric{distances};
    parlay::sequence<gbbs::uintE> centers = {1, 4, 6};
    auto state = gbbs::kmedian::BuildSolutionStateFromDistances(metric, centers);
    gbbs::kmedian::StaticSwapTable static_table(metric, state);
    gbbs::kmedian::DynamicSwapTable table(static_table, distances.size());

    auto apply = [&](gbbs::uintE p, gbbs::uintE q, double value) {
        const DistanceMatrix old = distances;
        auto old_distance = [&](gbbs::uintE u, gbbs::uintE v) { return old[u][v]; };
        distances[p][q] = distances[q][p] = value;
        auto new_distance = [&](gbbs::uintE u, gbbs::uintE v) { return distances[u][v]; };
        table.ApplyDistanceDecrease(p, q, state, old_distance, new_distance);
        VerifyWholeTable(table, ComputeBruteSwapTable(distances, centers, table));
    };

    apply(0, 7, 2.0);
    apply(1, 7, 1.0);
    apply(1, 6, 2.0);
    apply(4, 7, 0.5);
    apply(0, 5, 1.5);
    std::cout << "DynamicSwapTable sequential distance-decrease test passed." << std::endl;
}

void RunRandomizedDistanceDecreaseSequenceTest() {
    constexpr size_t kNumVertices = 10;
    constexpr size_t kNumAttempts = 500;
    constexpr uint64_t kSeed = 42;
    constexpr double kMinimumDistance = 1e-4;

    DistanceMatrix distances = MakePathDistances(kNumVertices);
    const MatrixMetric metric{distances};
    parlay::sequence<gbbs::uintE> centers = {1, 4, 7};
    auto state = gbbs::kmedian::BuildSolutionStateFromDistances(metric, centers);
    gbbs::kmedian::StaticSwapTable static_table(metric, state);
    gbbs::kmedian::DynamicSwapTable table(static_table, distances.size());

    std::mt19937_64 rng(kSeed);
    std::uniform_int_distribution<size_t> vertex_dist(0, kNumVertices - 1);
    std::uniform_real_distribution<double> factor_dist(0.05, 0.95);

    for (size_t attempt = 0; attempt < kNumAttempts; ++attempt) {
        gbbs::uintE p = static_cast<gbbs::uintE>(vertex_dist(rng));
        gbbs::uintE q = static_cast<gbbs::uintE>(vertex_dist(rng));
        while (q == p) {
            q = static_cast<gbbs::uintE>(vertex_dist(rng));
        }

        const double old_value = distances[p][q];
        if (old_value <= kMinimumDistance) {
            continue;
        }

        const double candidate = std::max(kMinimumDistance, old_value * factor_dist(rng));
        const double new_value = static_cast<double>(
            static_cast<gbbs::kmedian::Distance>(candidate));
        if (new_value >= old_value) {
            continue;
        }

        const DistanceMatrix old = distances;
        auto old_distance = [&](gbbs::uintE u, gbbs::uintE v) { return old[u][v]; };
        distances[p][q] = distances[q][p] = new_value;
        auto new_distance = [&](gbbs::uintE u, gbbs::uintE v) { return distances[u][v]; };
        table.ApplyDistanceDecrease(p, q, state, old_distance, new_distance);
        VerifyWholeTable(table, ComputeBruteSwapTable(distances, centers, table));
    }

    std::cout << "DynamicSwapTable randomized distance-decrease sequence test passed." << std::endl;
}

void RunRandomizedMutationTest() {
    const DistanceMatrix distances = MakePathDistances(6);
    const MatrixMetric metric{distances};
    parlay::sequence<gbbs::uintE> centers = {1, 4};
    auto state = gbbs::kmedian::BuildSolutionStateFromDistances(metric, centers);
    gbbs::kmedian::StaticSwapTable static_table(metric, state);
    gbbs::kmedian::DynamicSwapTable table(static_table, distances.size());

    std::vector<std::vector<double>> brute(
        static_table.num_rows(), std::vector<double>(static_table.num_cols()));
    for (size_t row = 0; row < static_table.num_rows(); ++row) {
        for (size_t col = 0; col < static_table.num_cols(); ++col) {
            brute[row][col] = static_table.gain(row, col);
        }
    }

    std::mt19937_64 rng(42);
    std::uniform_int_distribution<size_t> row_dist(0, table.num_rows() - 1);
    std::uniform_int_distribution<size_t> col_dist(0, table.num_cols() - 1);
    std::uniform_real_distribution<double> delta_dist(-10.0, 10.0);
    std::uniform_int_distribution<int> operation_dist(0, 1);

    for (size_t update = 0; update < 500; ++update) {
        const size_t col = col_dist(rng);
        const double delta = delta_dist(rng);
        if (operation_dist(rng) == 0) {
            const size_t row = row_dist(rng);
            table.PointAdd(col, row, delta);
            brute[row][col] += delta;
        } else {
            size_t left = row_dist(rng);
            size_t right = row_dist(rng);
            if (left > right) {
                std::swap(left, right);
            }
            table.RangeAdd(col, left, right, delta);
            for (size_t row = left; row <= right; ++row) {
                brute[row][col] += delta;
            }
        }
        VerifyWholeTable(table, brute);
    }

    std::cout << "DynamicSwapTable randomized mutation test passed." << std::endl;
}

}  // namespace

int main() {
    RunBasicDistanceDecreaseTests();
    RunOptimizedMixedDistanceDecreaseTests();
    RunSequentialDistanceDecreaseTest();
    RunRandomizedDistanceDecreaseSequenceTest();
    RunRandomizedMutationTest();
    return 0;
}
