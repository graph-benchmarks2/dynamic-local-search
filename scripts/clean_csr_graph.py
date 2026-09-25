from pathlib import Path
import sys

def clean(path_in, path_out):
    path_in = Path(path_in)
    path_out = Path(path_out)
    path_out.parent.mkdir(parents=True, exist_ok=True)

    print(f"Reading {path_in}")

    with path_in.open() as f:
        header = f.readline().strip()
        if header not in {"AdjacencyGraph", "WeightedAdjacencyGraph"}:
            raise ValueError(f"Unexpected header: {header}")

        weighted = header == "WeightedAdjacencyGraph"

        n = int(f.readline())
        m = int(f.readline())
        offsets = [int(f.readline()) for _ in range(n)]
        nbrs = [int(f.readline()) for _ in range(m)]
        weights = [int(f.readline()) for _ in range(m)] if weighted else None

    edges = {}

    raw_self_loops = 0
    raw_entries = 0
    bad_weights = 0

    for u in range(n):
        start = offsets[u]
        end = offsets[u + 1] if u + 1 < n else m

        for idx in range(start, end):
            v = nbrs[idx]
            raw_entries += 1

            if u == v:
                raw_self_loops += 1
                continue

            if v < 0 or v >= n:
                raise ValueError(f"Bad neighbor id {v} at edge index {idx}")

            a, b = (u, v) if u < v else (v, u)

            if weighted:
                w = weights[idx]
                if w <= 0:
                    bad_weights += 1
                old = edges.get((a, b))
                if old is None or w < old:
                    edges[(a, b)] = w
            else:
                edges[(a, b)] = None

    active = bytearray(n)
    for u, v in edges.keys():
        active[u] = 1
        active[v] = 1

    old_to_new = [-1] * n
    new_n = 0
    for u in range(n):
        if active[u]:
            old_to_new[u] = new_n
            new_n += 1

    if weighted:
        adj = [[] for _ in range(new_n)]
        for (u, v), w in edges.items():
            nu = old_to_new[u]
            nv = old_to_new[v]
            adj[nu].append((nv, w))
            adj[nv].append((nu, w))
    else:
        adj = [[] for _ in range(new_n)]
        for u, v in edges.keys():
            nu = old_to_new[u]
            nv = old_to_new[v]
            adj[nu].append(nv)
            adj[nv].append(nu)

    new_offsets = []
    new_m = 0

    if weighted:
        for u in range(new_n):
            adj[u].sort(key=lambda x: x[0])
            new_offsets.append(new_m)
            new_m += len(adj[u])
    else:
        for u in range(new_n):
            adj[u].sort()
            new_offsets.append(new_m)
            new_m += len(adj[u])

    print(f"Writing {path_out}")

    with path_out.open("w") as out:
        out.write(f"{header}\n")
        out.write(f"{new_n}\n")
        out.write(f"{new_m}\n")

        for off in new_offsets:
            out.write(f"{off}\n")

        if weighted:
            # GBBS weighted format: offsets, neighbors, weights.
            for u in range(new_n):
                for v, _ in adj[u]:
                    out.write(f"{v}\n")
            for u in range(new_n):
                for _, w in adj[u]:
                    out.write(f"{w}\n")
        else:
            for u in range(new_n):
                for v in adj[u]:
                    out.write(f"{v}\n")

    old_undirected_estimate = raw_entries // 2
    duplicate_or_asym_entries_removed = raw_entries - (2 * len(edges)) - raw_self_loops

    print(f"old_header={header}")
    print(f"old_n={n}")
    print(f"old_m={m}")
    print(f"new_n={new_n}")
    print(f"new_m={new_m}")
    print(f"distinct_undirected_edges={len(edges)}")
    print(f"removed_vertices={n - new_n}")
    print(f"removed_self_loop_entries={raw_self_loops}")
    print(f"removed_duplicate_or_asymmetric_entries={duplicate_or_asym_entries_removed}")
    if weighted:
        print(f"nonpositive_weight_entries_seen={bad_weights}")

    if new_m != 2 * len(edges):
        raise AssertionError("Output is not symmetric as expected.")
    if bad_weights:
        raise AssertionError("Weighted graph contains non-positive weights.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 clean_csr_graph.py <input.adj> <output.adj>")
        sys.exit(1)

    clean(sys.argv[1], sys.argv[2])
