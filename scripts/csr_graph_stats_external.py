from pathlib import Path
import argparse
import subprocess
import tempfile
import os
import shutil

def run(cmd):
    return subprocess.check_output(cmd, text=True).strip()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("graph")
    parser.add_argument("--tmpdir", default=None)
    args = parser.parse_args()

    path = Path(args.graph)

    tmpdir = Path(args.tmpdir) if args.tmpdir else Path(tempfile.mkdtemp(prefix="csr_stats_"))
    tmpdir.mkdir(parents=True, exist_ok=True)

    canonical_edges = tmpdir / "canonical_edges.txt"

    try:
        with path.open() as f, canonical_edges.open("w") as out:
            header = f.readline().strip()
            weighted = header == "WeightedAdjacencyGraph"

            if header not in {"AdjacencyGraph", "WeightedAdjacencyGraph"}:
                raise RuntimeError(f"Unexpected header: {header}")

            n = int(f.readline())
            m = int(f.readline())

            offsets = [int(f.readline()) for _ in range(n)]

            isolated_vertices = 0
            self_loop_entries = 0
            non_loop_adjacency_entries = 0

            for u in range(n):
                start = offsets[u]
                end = offsets[u + 1] if u + 1 < n else m
                deg = end - start

                if deg == 0:
                    isolated_vertices += 1

                for _ in range(deg):
                    v = int(f.readline())

                    if u == v:
                        self_loop_entries += 1
                    else:
                        non_loop_adjacency_entries += 1
                        a, b = (u, v) if u < v else (v, u)
                        out.write(f"{a} {b}\n")

        sorted_unique = tmpdir / "canonical_edges_unique.txt"

        subprocess.check_call(
            f"LC_ALL=C sort -S 50% -T {tmpdir} -u {canonical_edges} > {sorted_unique}",
            shell=True,
        )

        distinct_undirected_edges = int(run(["wc", "-l", str(sorted_unique)]).split()[0])

        # non_loop_adjacency_entries == 2 * distinct_undirected_edges
        excess_or_missing = non_loop_adjacency_entries - 2 * distinct_undirected_edges

        print(f"file: {path}")
        print(f"header: {header}")
        print(f"nodes_total: {n}")
        print(f"adjacency_entries_total: {m}")
        print(f"isolated_vertices: {isolated_vertices}")
        print(f"self_loop_entries: {self_loop_entries}")
        print(f"non_loop_adjacency_entries: {non_loop_adjacency_entries}")
        print(f"distinct_undirected_edges_no_self_loops: {distinct_undirected_edges}")
        print(f"expected_clean_symmetric_adjacency_entries: {2 * distinct_undirected_edges}")
        print(f"duplicate_or_asymmetry_excess_entries: {excess_or_missing}")

        if self_loop_entries == 0 and isolated_vertices == 0 and excess_or_missing == 0:
            print("status: CLEAN_SIMPLE_SYMMETRIC_NO_ISOLATES")
        else:
            print("status: NEEDS_CLEANING_OR_CHECK")

    finally:
        if args.tmpdir is None:
            shutil.rmtree(tmpdir, ignore_errors=True)

if __name__ == "__main__":
    main()
