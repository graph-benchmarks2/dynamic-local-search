#!/usr/bin/env python3

import argparse
import random
from pathlib import Path


def canonical_edge(u: int, v: int) -> tuple[int, int]:
    return (u, v) if u < v else (v, u)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a connected undirected unweighted graph plus a sequence "
            "of distinct edge insertions that are absent from the initial graph."
        )
    )
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--initial-m", type=int, default=250)
    parser.add_argument("--updates", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--graph-out",
        type=Path,
        default=Path("inputs/DynamicStress/stress_n100_m250.edgelist"),
    )
    parser.add_argument(
        "--updates-out",
        type=Path,
        default=Path("inputs/DynamicStress/stress_n100_updates100.txt"),
    )
    args = parser.parse_args()

    if args.n < 2:
        raise SystemExit("ERROR: --n must be at least 2.")

    max_edges = args.n * (args.n - 1) // 2

    if args.initial_m < args.n - 1:
        raise SystemExit(
            "ERROR: --initial-m must be at least n-1 so the graph can be connected."
        )

    if args.initial_m + args.updates > max_edges:
        raise SystemExit(
            "ERROR: initial edges plus updates exceed the number of possible edges."
        )

    rng = random.Random(args.seed)
    edges: set[tuple[int, int]] = set()

    # Random spanning tree: every new vertex attaches to one earlier vertex.
    # This guarantees connectedness without forcing a path-shaped graph.
    order = list(range(args.n))
    rng.shuffle(order)

    for i in range(1, args.n):
        u = order[i]
        v = order[rng.randrange(i)]
        edges.add(canonical_edge(u, v))

    # Fill the initial graph with additional random distinct edges.
    while len(edges) < args.initial_m:
        u = rng.randrange(args.n)
        v = rng.randrange(args.n)

        if u == v:
            continue

        edges.add(canonical_edge(u, v))

    initial_edges = sorted(edges)

    # Every update is a distinct edge that is absent at the time it is chosen.
    updates: list[tuple[int, int]] = []

    while len(updates) < args.updates:
        u = rng.randrange(args.n)
        v = rng.randrange(args.n)

        if u == v:
            continue

        edge = canonical_edge(u, v)

        if edge in edges:
            continue

        edges.add(edge)
        updates.append(edge)

    args.graph_out.parent.mkdir(parents=True, exist_ok=True)
    args.updates_out.parent.mkdir(parents=True, exist_ok=True)

    with args.graph_out.open("w", encoding="utf-8") as out:
        for u, v in initial_edges:
            out.write(f"{u} {v}\n")

    with args.updates_out.open("w", encoding="utf-8") as out:
        for u, v in updates:
            out.write(f"{u} {v}\n")

    print(f"n = {args.n}")
    print(f"initial_m = {len(initial_edges)}")
    print(f"updates = {len(updates)}")
    print(f"graph = {args.graph_out}")
    print(f"update_sequence = {args.updates_out}")


if __name__ == "__main__":
    main()

