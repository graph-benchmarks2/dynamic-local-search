#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>

#include "parlay/sequence.h"

namespace gbbs {
namespace kmedian {

class RangeMaxTree {
public:
    RangeMaxTree() = default;

    explicit RangeMaxTree(const parlay::sequence<double>& values) {
        Build(values);
    }

    RangeMaxTree(
        const parlay::sequence<double>& values,
        const size_t* tie_break_keys) {

        Build(values, tie_break_keys);
    }

    void Build(const parlay::sequence<double>& values) {
        Build(values, nullptr);
    }

    void Build(
        const parlay::sequence<double>& values,
        const size_t* tie_break_keys) {

        size_ = values.size();
        tie_break_keys_ = tie_break_keys;

        if (size_ == 0) {
            std::cout << "ERROR: RangeMaxTree requires at least one value." << std::endl;
            std::exit(-1);
        }

        max_value_ = parlay::sequence<double>(
            4 * size_,
            -std::numeric_limits<double>::infinity());

        max_index_ = parlay::sequence<size_t>(
            4 * size_,
            std::numeric_limits<size_t>::max());

        lazy_add_ = parlay::sequence<double>(4 * size_, 0.0);

        BuildRecursive(1, 0, size_ - 1, values);
    }

    size_t size() const {
        return size_;
    }

    void RangeAdd(size_t left, size_t right, double delta) {
        if (left > right) {
            return;
        }

        CheckIndex(left);
        CheckIndex(right);

        RangeAddRecursive(1, 0, size_ - 1, left, right, delta);
    }

    void PointAdd(size_t index, double delta) {
        CheckIndex(index);
        RangeAdd(index, index, delta);
    }

    double MaxValue() const {
        CheckBuilt();
        return max_value_[1];
    }

    size_t MaxIndex() const {
        CheckBuilt();
        return max_index_[1];
    }

    double Value(size_t index) {
        CheckIndex(index);
        return ValueRecursive(1, 0, size_ - 1, index);
    }

private:
    void CheckBuilt() const {
        if (size_ == 0) {
            std::cout << "ERROR: RangeMaxTree has not been initialized." << std::endl;
            std::exit(-1);
        }
    }

    void CheckIndex(size_t index) const {
        CheckBuilt();

        if (index >= size_) {
            std::cout << "ERROR: RangeMaxTree index " << index
                      << " is outside range [0, " << size_ << ")."
                      << std::endl;
            std::exit(-1);
        }
    }

    bool Better(
        double candidate_value,
        size_t candidate_index,
        double current_value,
        size_t current_index) const {

        if (candidate_value > current_value) {
            return true;
        }

        if (candidate_value < current_value) {
            return false;
        }

        if (tie_break_keys_ != nullptr) {
            return tie_break_keys_[candidate_index]
                < tie_break_keys_[current_index];
        }

        return candidate_index < current_index;
    }

    void Apply(size_t node, double delta) {
        max_value_[node] += delta;
        lazy_add_[node] += delta;
    }

    void Push(size_t node) {
        if (lazy_add_[node] == 0.0) {
            return;
        }

        Apply(node * 2, lazy_add_[node]);
        Apply(node * 2 + 1, lazy_add_[node]);

        lazy_add_[node] = 0.0;
    }

    void Pull(size_t node) {
        const size_t left_child = node * 2;
        const size_t right_child = node * 2 + 1;

        if (Better(
                max_value_[left_child],
                max_index_[left_child],
                max_value_[right_child],
                max_index_[right_child])) {

            max_value_[node] = max_value_[left_child];
            max_index_[node] = max_index_[left_child];
        } else {
            max_value_[node] = max_value_[right_child];
            max_index_[node] = max_index_[right_child];
        }
    }

    void BuildRecursive(
        size_t node,
        size_t left,
        size_t right,
        const parlay::sequence<double>& values) {

        if (left == right) {
            max_value_[node] = values[left];
            max_index_[node] = left;
            return;
        }

        const size_t middle = left + (right - left) / 2;

        BuildRecursive(node * 2, left, middle, values);
        BuildRecursive(node * 2 + 1, middle + 1, right, values);

        Pull(node);
    }

    void RangeAddRecursive(
        size_t node,
        size_t left,
        size_t right,
        size_t query_left,
        size_t query_right,
        double delta) {

        if (query_left <= left && right <= query_right) {
            Apply(node, delta);
            return;
        }

        Push(node);

        const size_t middle = left + (right - left) / 2;

        if (query_left <= middle) {
            RangeAddRecursive(
                node * 2,
                left,
                middle,
                query_left,
                query_right,
                delta);
        }

        if (query_right > middle) {
            RangeAddRecursive(
                node * 2 + 1,
                middle + 1,
                right,
                query_left,
                query_right,
                delta);
        }

        Pull(node);
    }

    double ValueRecursive(
        size_t node,
        size_t left,
        size_t right,
        size_t index) {

        if (left == right) {
            return max_value_[node];
        }

        Push(node);

        const size_t middle = left + (right - left) / 2;

        if (index <= middle) {
            return ValueRecursive(node * 2, left, middle, index);
        }

        return ValueRecursive(node * 2 + 1, middle + 1, right, index);
    }

    size_t size_ = 0;
    const size_t* tie_break_keys_ = nullptr;

    parlay::sequence<double> max_value_;
    parlay::sequence<size_t> max_index_;
    parlay::sequence<double> lazy_add_;
};

}  // namespace kmedian
}  // namespace gbbs
