#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

#include "benchmarks/Clustering/K-Median/common/RangeMaxTree.h"
#include "parlay/sequence.h"

namespace {

void CheckEqual(double actual, double expected, const char* label) {
    constexpr double kTolerance = 1e-12;

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

void RunRandomizedStressTest() {
    constexpr size_t kSize = 32;
    constexpr size_t kNumUpdates = 5000;
    constexpr uint64_t kSeed = 42;

    std::mt19937_64 rng(kSeed);

    std::uniform_real_distribution<double> value_dist(-100.0, 100.0);
    std::uniform_real_distribution<double> delta_dist(-20.0, 20.0);
    std::uniform_int_distribution<size_t> index_dist(0, kSize - 1);

    parlay::sequence<double> initial_values(kSize);
    std::vector<double> brute_values(kSize);

    for (size_t i = 0; i < kSize; ++i) {
        const double value = value_dist(rng);
        initial_values[i] = value;
        brute_values[i] = value;
    }

    gbbs::kmedian::RangeMaxTree tree(initial_values);

    for (size_t update = 0; update < kNumUpdates; ++update) {
        size_t left = index_dist(rng);
        size_t right = index_dist(rng);

        if (left > right) {
            std::swap(left, right);
        }

        const double delta = delta_dist(rng);

        tree.RangeAdd(left, right, delta);

        for (size_t i = left; i <= right; ++i) {
            brute_values[i] += delta;
        }

        for (size_t i = 0; i < kSize; ++i) {
            CheckEqual(
                tree.Value(i),
                brute_values[i],
                "randomized point value");
        }

        double expected_max = brute_values[0];
        size_t expected_index = 0;

        for (size_t i = 1; i < kSize; ++i) {
            if (brute_values[i] > expected_max ||
                (brute_values[i] == expected_max && i < expected_index)) {

                expected_max = brute_values[i];
                expected_index = i;
            }
        }

        CheckEqual(
            tree.MaxValue(),
            expected_max,
            "randomized max value");

        CheckEqual(
            tree.MaxIndex(),
            expected_index,
            "randomized max index");
    }

    std::cout << "RangeMaxTree randomized stress test passed." << std::endl;
}

}  // namespace

int main() {
    parlay::sequence<double> values = {3.0, -1.0, 7.0, 2.0, 5.0};

    gbbs::kmedian::RangeMaxTree tree(values);

    CheckEqual(tree.MaxValue(), 7.0, "initial max value");
    CheckEqual(tree.MaxIndex(), size_t{2}, "initial max index");

    for (size_t i = 0; i < values.size(); ++i) {
        CheckEqual(tree.Value(i), values[i], "initial point value");
    }

    tree.RangeAdd(1, 3, 4.0);

    // [3, 3, 11, 6, 5]
    CheckEqual(tree.MaxValue(), 11.0, "max after first range add");
    CheckEqual(tree.MaxIndex(), size_t{2}, "max index after first range add");
    CheckEqual(tree.Value(1), 3.0, "value 1 after first range add");
    CheckEqual(tree.Value(2), 11.0, "value 2 after first range add");
    CheckEqual(tree.Value(3), 6.0, "value 3 after first range add");

    tree.PointAdd(4, 10.0);

    // [3, 3, 11, 6, 15]
    CheckEqual(tree.MaxValue(), 15.0, "max after point add");
    CheckEqual(tree.MaxIndex(), size_t{4}, "max index after point add");
    CheckEqual(tree.Value(4), 15.0, "value 4 after point add");

    tree.RangeAdd(0, 4, -5.0);

    // [-2, -2, 6, 1, 10]
    CheckEqual(tree.MaxValue(), 10.0, "max after full range add");
    CheckEqual(tree.MaxIndex(), size_t{4}, "max index after full range add");

    tree.PointAdd(0, 12.0);

    // [10, -2, 6, 1, 10]
    //
    // Tie between indices 0 and 4. We deliberately choose the smaller index.
    CheckEqual(tree.MaxValue(), 10.0, "max after tie");
    CheckEqual(tree.MaxIndex(), size_t{0}, "tie-breaking index");

    tree.RangeAdd(1, 4, 4.0);

    // [10, 2, 10, 5, 14]
    CheckEqual(tree.MaxValue(), 14.0, "final max value");
    CheckEqual(tree.MaxIndex(), size_t{4}, "final max index");

    const parlay::sequence<double> expected = {10.0, 2.0, 10.0, 5.0, 14.0};

    for (size_t i = 0; i < expected.size(); ++i) {
        CheckEqual(tree.Value(i), expected[i], "final point value");
    }

    std::cout << "RangeMaxTree test passed." << std::endl;

    RunRandomizedStressTest();

    return 0;
}
