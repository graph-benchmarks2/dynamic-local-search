#!/usr/bin/env python3
from pathlib import Path
import argparse


MASK64 = (1 << 64) - 1


def splitmix64(x: int) -> int:
    x = (x + 0x9E3779B97F4A7C15) & MASK64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & MASK64
    return (x ^ (x >> 31)) & MASK64


def edge_weight(u: int, v: int, seed: int, max_weight: int) -> int:
    a, b = (u, v) if u <= v else (v, u)

    # Combine endpoints and seed deterministically.
    x = seed & MASK64
    x ^= (a + 0x9E3779B97F4A7C15) & MASK64
    x = splitmix64(x)
    x ^= (b + 0xBF58476D1CE4E5B9) & MASK64
    x = splitmix64(x)

    return 1 + (x % max_weight)


def read_int_line(f, what: str) -> int:
    line = f.readline()
    if line == "":
        raise RuntimeError(f"Unexpected EOF while reading {what}")
    return int(line.strip())


def main():
    parser = argparse.ArgumentParser(
        description="Convert an unweighted GBBS CSR .adj graph to weighted CSR "
                    "with deterministic uniform integer weights."
    )
    parser.add_argument("input", help="Input unweighted GBBS .adj file")
    parser.add_argument("output", help="Output weighted GBBS .adj file")
    parser.add_argument("--max-weight", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--force", action="store_true",
                        help="Allow overwriting existing output")
    args = parser.parse_args()

    path_in = Path(args.input)
    path_out = Path(args.output)

    if args.max_weight <= 0:
        raise ValueError("--max-weight must be positive")

    if path_out.exists() and not args.force:
        raise FileExistsError(f"Output exists: {path_out}. Use --force to overwrite.")

    path_out.parent.mkdir(parents=True, exist_ok=True)

    with path_in.open("r") as f, path_out.open("w") as out:
        header = f.readline().strip()

        if header != "AdjacencyGraph":
            raise RuntimeError(
                f"Expected unweighted AdjacencyGraph input, got: {header}"
            )

        n = read_int_line(f, "n")
        m = read_int_line(f, "m")

        offsets = [read_int_line(f, f"offset[{i}]") for i in range(n)]

        out.write("WeightedAdjacencyGraph\n")
        out.write(f"{n}\n")
        out.write(f"{m}\n")

        for off in offsets:
            out.write(f"{off}\n")

        # Stream neighbors to output and store them temporarily.
        # For huge graphs this list is large, but comparable to the input's
        # adjacency block and avoids rereading the file.
        nbrs = [0] * m

        for i in range(m):
            v = read_int_line(f, f"neighbor[{i}]")
            if v < 0 or v >= n:
                raise RuntimeError(f"Invalid neighbor id {v} at adjacency index {i}")
            nbrs[i] = v
            out.write(f"{v}\n")

        # Append weight block in the exact same adjacency-entry order.
        for u in range(n):
            start = offsets[u]
            end = offsets[u + 1] if u + 1 < n else m

            for idx in range(start, end):
                v = nbrs[idx]
                w = edge_weight(u, v, args.seed, args.max_weight)
                out.write(f"{w}\n")

    print(f"input: {path_in}")
    print(f"output: {path_out}")
    print(f"nodes: {n}")
    print(f"adjacency_entries: {m}")
    print(f"max_weight: {args.max_weight}")
    print(f"seed: {args.seed}")
    print("status: OK")


if __name__ == "__main__":
    main()
