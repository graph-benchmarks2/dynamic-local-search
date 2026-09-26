# Dynamic Local Search for Graph k-Median

This repository contains the C++20 implementations used for the experimental evaluation of distance-decrease local search for the $k$-median objective. The implementations are integrated into the [Graph-Based Benchmark Suite (GBBS)](https://github.com/ParAlg/gbbs).

## 1. Implemented Algorithms

The evaluation considers three local search strategies:

* **FastPAM1-Cold:** Computes a locally optimal solution from scratch using a random initial set of $k$ centers. At every iteration, it performs the best improving swap in the entire swap table.
* **FastPAM1-Warm:** Starts from the center set obtained for the previous graph state but reconstructs its auxiliary data structures, including the swap table, after each edge insertion.
* **Dynamic $\mathcal{LS}$:** Maintains both the center set and its dynamic swap table across edge insertions. The dynamic range data structures are segment trees with lazy propagation.

All strategies accept every strictly improving swap; no swap-gain threshold is applied.

## 2. Building

The implementations are part of GBBS and use C++20. The experiments described in the paper were compiled using GCC 12 through Bazel with `-O2` and `-march=native`.

From the repository root, build the three executables:

```bash
bazel build \
  //benchmarks/Clustering/K-Median/StaticLocalSearch:StaticLocalSearch_main \
  //benchmarks/Clustering/K-Median/StaticLocalSearch:WarmStartStaticLocalSearch_main \
  //benchmarks/Clustering/K-Median/DynamicLocalSearch:DynamicLocalSearch_main
```

The resulting executables are located in the corresponding directories under `bazel-bin/`.

## 3. Datasets

We use two temporal networks from the [higher-order network collection of Benson et al.](https://www.cs.cornell.edu/~arb/data/):

* **DBLP:** A co-authorship network in which vertices represent authors and edges represent collaborations.
* **Stack Overflow:** A co-occurrence network in which vertices represent tags and edges connect tags appearing in the same thread.

The source datasets contain timestamped simplices. We project each simplex onto all unordered pairs of its vertices to obtain an undirected, unweighted temporal graph. Only the first occurrence of each edge is retained, such that repeated co-occurrences do not generate additional insertions.

### 3.1. Preparing the temporal snapshots

For each dataset and selected cutoff time $T$, the initial graph and its update sequence are constructed as follows:

1. Process the temporal data in chronological order and construct the graph containing all distinct edges observed up to and including $T$.
2. Compute the largest connected component (LCC) of this graph and retain only its vertices and induced edges as the initial graph.
3. Continue processing the temporal data after $T$, retaining an edge insertion only when both endpoints belong to the fixed LCC and the edge has not previously occurred.
4. Preserve the chronological order of the retained insertions.

Thus, each instance has a fixed vertex set and an update sequence consisting exclusively of genuine edge insertions. No multi-edges or new vertices are introduced at any point in time.

We use DBLP snapshots at the years 1977, 1983, 1984 and 1990. For Stack Overflow, snapshots are taken after approximately 5%, 40% and 80% of its temporal evolution. For each snapshot, we consider prefixes of the resulting update sequence containing 100, 500 and 1,000 insertions.

### 3.2. Initial graph sizes

| Dataset        | Cutoff | Vertices $n$ | Edges $m$ | $m/n$ |
| -------------- | -----: | -----------: | --------: | ----: |
| DBLP           |   1977 |        5,263 |    10,319 |  1.96 |
| DBLP           |   1983 |       16,134 |    33,413 |  2.07 |
| DBLP           |   1984 |       19,432 |    40,425 |  2.08 |
| DBLP           |   1990 |       55,067 |   122,282 |  2.22 |
| Stack Overflow |     5% |       22,126 |   520,820 | 23.54 |
| Stack Overflow |    40% |       35,397 | 2,179,797 | 61.58 |
| Stack Overflow |    80% |       45,497 | 3,505,795 | 77.06 |

The original datasets are available from:

* [DBLP co-authorship data](https://www.cs.cornell.edu/~arb/data/coauth-DBLP/)
* [Stack Overflow tag data](https://www.cs.cornell.edu/~arb/data/tags-stack-overflow/)

The original datasets and generated graph instances are not included in this repository.

## 4. Experimental Configuration

We consider $k \in {10,20,30,50,100}$ centers and update sequences of 100, 500 and 1,000 edge insertions. FastPAM1-Cold is evaluated at the corresponding checkpoints while FastPAM1-Warm and Dynamic $\mathcal{LS}$ process insertions incrementally.

The experiments described in the paper were run on Debian 12 using a virtual machine on a server equipped with two AMD EPYC 7763 64-core processors and approximately 4 TB of available RAM. Experiments were repeated three times and median running times were reported. The timeout for each individual experiment was 2.5 hours per run.

All strategies use exact all-pairs shortest-path (APSP) distances. Shortest-path computations and APSP repairs are excluded from the reported local search running times.

## 5. Executing the Experiments

The following examples use $k=10$, random initialization, seed 42 and `-ls-delta 0`.

Replace the example input paths with the locations of your prepared graph and update files.

**FastPAM1-Cold**

```bash
bazel-bin/benchmarks/Clustering/K-Median/StaticLocalSearch/StaticLocalSearch_main \
  -k 10 -seed 42 -init random -ls-delta 0 \
  inputs/initial.adj
```

**FastPAM1-Warm**

```bash
bazel-bin/benchmarks/Clustering/K-Median/StaticLocalSearch/WarmStartStaticLocalSearch_main \
  -k 10 -seed 42 -init random -ls-delta 0 \
  -updates inputs/updates_100.txt \
  inputs/initial.adj
```

**Dynamic $\mathcal{LS}$**

```bash
bazel-bin/benchmarks/Clustering/K-Median/DynamicLocalSearch/DynamicLocalSearch_main \
  -k 10 -seed 42 -init random -ls-delta 0 \
  -updates inputs/updates_100.txt \
  inputs/initial.adj
```

The graph files use a format supported by GBBS. For the unweighted graphs in our experiments, the update file is a whitespace-separated list of vertex pairs with one insertion per line:

```text
12 37
8 19
37 42
```

Vertex IDs must be valid for the initial graph. The update file must contain genuine edge insertions in chronological order without self-loops or duplicate edges.

The incremental implementations initialize a solution on the initial graph and then process the supplied update sequence. 

### 5.1. Optional Cold-Start Quality Comparisons

To reproduce the FastPAM1-Cold checkpoint comparisons there is two options: 

1) Apply the corresponding update prefix to the initial graph and run the static executable on the resulting graph. This is particularly useful for running time comparisons.
2) Both FastPAM1-Warm and Dynamic $\mathcal{LS}$ support optional cold-start quality evaluations at selected update checkpoints: 

These run FastPAM1-Cold on the graph after the specified number of insertions, allowing direct comparison of the resulting $k$-median objective values.

For example, the following command runs Dynamic $\mathcal{LS}$ on 1,000 insertions and additionally evaluates cold-start solutions at checkpoints 100, 500 and 1,000:

```bash
bazel-bin/benchmarks/Clustering/K-Median/DynamicLocalSearch/DynamicLocalSearch_main \
  -k 10 -seed 42 -init random -ls-delta 0 \
  -updates inputs/updates_1000.txt \
  -quality-cold-trials 3 \
  -quality-cold-seed 42 \
  -quality-cold-checkpoints 100,500,1000 \
  inputs/initial.adj
```

The same quality-evaluation flags are supported by `WarmStartStaticLocalSearch_main`. Set `-quality-cold-trials 0` (the default) to disable these additional evaluations. The cold-start quality comparisons are separate from the reported incremental local-search update times.

## 6. Acknowledgments

This implementation builds on the Graph-Based Benchmark Suite (GBBS). The original GBBS license and attribution are preserved in this repository.

The temporal datasets are provided by Benson et al., *Simplicial Closure and Higher-Order Link Prediction*, PNAS, 2018.

