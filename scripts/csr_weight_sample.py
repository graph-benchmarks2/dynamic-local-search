#!/usr/bin/env python3
from pathlib import Path
import argparse
import random
import statistics


def read_int_line(f, what):
    line = f.readline()
    if line == "":
        raise RuntimeError(f"Unexpected EOF while reading {what}")
    return int(line.strip())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("graph")
    parser.add_argument("--sample", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--hist", action="store_true")
    args = parser.parse_args()

    path = Path(args.graph)
    rng = random.Random(args.seed)

    with path.open() as f:
        header = f.readline().strip()

        if header not in {"AdjacencyGraph", "WeightedAdjacencyGraph"}:
            raise RuntimeError(f"Unexpected header: {header}")

        weighted = header == "WeightedAdjacencyGraph"

        n = read_int_line(f, "n")
        m = read_int_line(f, "m")

        offsets = [read_int_line(f, f"offset[{i}]") for i in range(n)]

        # Skip adjacency entries.
        for i in range(m):
            _ = f.readline()
            if _ == "":
                raise RuntimeError(f"Unexpected EOF while skipping adjacency entry {i}")

        print(f"file: {path}")
        print(f"header: {header}")
        print(f"weighted: {weighted}")
        print(f"nodes: {n}")
        print(f"adjacency_entries: {m}")

        if not weighted:
            print("status: unweighted graph, no weight block present")
            return

        sample_indices = set(rng.sample(range(m), min(args.sample, m)))
        sampled_weights = []
        min_w = None
        max_w = None
        sum_w = 0
        count = 0
        small_hist = {}

        for i in range(m):
            w = read_int_line(f, f"weight[{i}]")

            if i in sample_indices:
                sampled_weights.append(w)

            min_w = w if min_w is None else min(min_w, w)
            max_w = w if max_w is None else max(max_w, w)
            sum_w += w
            count += 1

            if args.hist and w <= 1000:
                small_hist[w] = small_hist.get(w, 0) + 1

        avg_w = sum_w / count if count else 0

        print(f"weight_entries: {count}")
        print(f"weight_min: {min_w}")
        print(f"weight_max: {max_w}")
        print(f"weight_avg: {avg_w:.6f}")

        if sampled_weights:
            sampled_sorted = sorted(sampled_weights)
            print(f"sample_count: {len(sampled_weights)}")
            print(f"sample_min: {min(sampled_weights)}")
            print(f"sample_max: {max(sampled_weights)}")
            print(f"sample_avg: {statistics.mean(sampled_weights):.6f}")
            print("sample_weights_sorted:")
            print(" ".join(map(str, sampled_sorted)))

        if args.hist:
            print("histogram_weights_leq_1000:")
            for w in sorted(small_hist):
                print(f"{w}\t{small_hist[w]}")


if __name__ == "__main__":
    main()
