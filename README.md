# GBBS: Graph Based Benchmark Suite  ![Bazel build](https://github.com/paralg/gbbs/workflows/CI/badge.svg)

Organization
--------

This repository contains code for our SPAA paper "Theoretically Efficient
Parallel Graph Algorithms Can Be Fast and Scalable" (SPAA'18). It includes
implementations of the following parallel graph algorithms:

**Clustering Problems**
* SCAN Graph Clustering
* Graph-Based Hierarchical Agglomerative Clustering (Graph HAC)

**Connectivity Problems**
* Low-Diameter Decomposition
* Connectivity
* Spanning Forest
* Biconnectivity
* Minimum Spanning Tree
* Strongly Connected Components

**Covering Problems**
* Coloring
* Maximal Matching
* Maximal Independent Set
* Approximate Set Cover

**Eigenvector Problems**
* PageRank

**Substructure Problems**
* Triangle Counting
* Approximate Densest Subgraph
* k-Core (Coreness)
* Degeneracy Ordering (Low-Outdegree Orientation)
* k-Clique Counting
* 5-Cycle Counting
* k-Truss

**Shortest Path Problems**
* Unweighted SSSP (Breadth-First Search)
* General Weight SSSP (Bellman-Ford)
* Integer Weight SSSP (Weighted Breadth-First Search)
* Single-Source Betweenness Centrality
* Single-Source Widest Path
* k-Spanner

The code for these applications is located in the `benchmark` directory. The
implementations are based on the Ligra/Ligra+/Julienne graph processing
frameworks. The framework code is located in the `src` directory.

If you use our work, please cite our [paper](https://arxiv.org/abs/1805.05208):

```
@inproceedings{dhulipala2018theoretically,
  author    = {Laxman Dhulipala and
               Guy E. Blelloch and
               Julian Shun},
  title     = {Theoretically Efficient Parallel Graph Algorithms Can Be Fast and
               Scalable},
  booktitle = {ACM Symposium on Parallelism in Algorithms and Architectures (SPAA)},
  year      = {2018},
}
```

Compilation
--------

Compiler:
* g++ &gt;= 7.4.0 with support for Cilk Plus
* g++ &gt;= 7.4.0 with pthread support (Homemade Scheduler)

Build system:
* [Bazel](https://docs.bazel.build/versions/master/install.html) 2.1.0
* Make --- though our primary build system is Bazel, we also maintain Makefiles
  for those who wish to run benchmarks without installing Bazel.

The default compilation uses a lightweight scheduler developed at CMU (Homemade)
for parallelism, which results in comparable performance to Cilk Plus. The
half-lengths for certain functions such as histogramming are lower using
Homemade, which results in better performance for codes like KCore.

The benchmark supports both uncompressed and compressed graphs. The uncompressed
format is identical to the uncompressed format in Ligra. The compressed format,
called bytepd_amortized (bytepda) is similar to the parallelByte format used in
Ligra+, with some additional functionality to support efficiently packs,
filters, and other operations over neighbor lists.

To compile codes for graphs with more than 2^32 edges, the `GBBSLONG` command-line
parameter should be set. If the graph has more than 2^32 vertices, the
`GBBSEDGELONG` command-line parameter should be set. Note that the codes have not
been tested with more than 2^32 vertices, so if any issues arise please contact
[Laxman Dhulipala](mailto:laxman@umd.edu).

To compile with the Cilk Plus scheduler instead of the Homegrown scheduler, use
the Bazel configuration `--config=cilk`. To compile using OpenMP instead, use
the Bazel configuration `--config=openmp`. To compile serially instead, use the
Bazel configuration `--config=serial`. (For the Makefiles, instead set the
environment variables `CILK`, `OPENMP`, or `SERIAL` respectively.)

To build:
```sh
# Load external libraries as submodules. (This only needs to be run once.)
git submodule update --init

# For Bazel:
$ bazel build  //...  # compiles all benchmarks

# For Make:
# First set the appropriate environment variables, e.g., first run
# `export CILK=1` to compile with Cilk Plus.
# After that, build using `make`.
$ cd benchmarks/BFS/NonDeterministicBFS  # go to a benchmark
$ make
```
Note that the default compilation mode in bazel is to build optimized binaries
(stripped of debug symbols). You can compile debug binaries by supplying `-c
dbg` to the bazel build command.

The following commands cleans the directory:
```sh
# For Bazel:
$ bazel clean  # removes all executables

# For Make:
$ make clean  # removes executables for the current directory
```

Running code
-------
The applications take the input graph as input as well as an optional
flag "-s" to indicate a symmetric graph.  Symmetric graphs should be
called with the "-s" flag for better performance. For example:

```sh
# For Bazel:
$ bazel run //benchmarks/BFS/NonDeterministicBFS:BFS_main -- -s -src 10 ~/gbbs/inputs/rMatGraph_J_5_100
$ bazel run //benchmarks/IntegralWeightSSSP/JulienneDBS17:wBFS_main -- -s -w -src 15 ~/gbbs/inputs/rMatGraph_WJ_5_100

# For Make:
$ ./BFS -s -src 10 ../../../inputs/rMatGraph_J_5_100
$ ./wBFS -s -w -src 15 ../../../inputs/rMatGraph_WJ_5_100
```

Note that the codes that compute single-source shortest paths (or centrality)
take an extra `-src` flag. The benchmark is run four times by default, and can
be changed by passing the `-rounds` flag followed by an integer indicating the
number of runs.

On NUMA machines, adding the command "numactl -i all " when running
the program may improve performance for large graphs. For example:

```sh
$ numactl -i all bazel run [...]
```

Running code on compressed graphs
-----------

We make use of the bytePDA format in our benchmark, which is similar to the
parallelByte format of Ligra+, extended with additional functionality. We have
provided a converter utility which takes as input an uncompressed graph and
outputs a bytePDA graph. The converter can be used as follows:

```sh
# For Bazel:
bazel run //utils:compressor -- -s -o ~/gbbs/inputs/rMatGraph_J_5_100.bytepda ~/gbbs/inputs/rMatGraph_J_5_100
bazel run //utils:compressor -- -s -w -o ~/gbbs/inputs/rMatGraph_WJ_5_100.bytepda ~/gbbs/inputs/rMatGraph_WJ_5_100

# For Make:
./compressor -s -o ../inputs/rMatGraph_J_5_100.bytepda ../inputs/rMatGraph_J_5_100
./compressor -s -w -o ../inputs/rMatGraph_WJ_5_100.bytepda ../inputs/rMatGraph_WJ_5_100
```

After an uncompressed graph has been converted to the bytepda format,
applications can be run on it by passing in the usual command-line flags, with
an additional `-c` flag.

```sh
# For Bazel:
$ bazel run //benchmarks/BFS/NonDeterministicBFS:BFS_main -- -s -c -src 10 ~/gbbs/inputs/rMatGraph_J_5_100.bytepda

# For Make:
$ ./BFS -s -c -src 10 ../../../inputs/rMatGraph_J_5_100.bytepda
$ ./wBFS -s -w -c -src 15 ../../../inputs/rMatGraph_WJ_5_100.bytepda
```

When processing large compressed graphs, using the `-m` command-line flag can
help if the file is already in the page cache, since the compressed graph data
can be mmap'd. Application performance will be affected if the file is not
already in the page-cache. We have found that using `-m` when the compressed
graph is backed by SSD results in a slow first-run, followed by fast subsequent
runs.

Running code on binary-encoded graphs
-----------
We make use of a binary-graph format in our benchmark. The binary representation
stores the representation we use for in-memory processing (compressed sparse row)
directly on disk, which enables applications to avoid string-conversion overheads
associated with the adjacency graph format described below. We have provided a
converter utility which takes as input an uncompressed graph (e.g., in adjacency
graph format) and outputs this graph in the binary format. The converter can be
used as follows:

```sh
# For Bazel:
bazel run //utils:compressor -- -s -o ~/gbbs/inputs/rMatGraph_J_5_100.binary ~/gbbs/inputs/rMatGraph_J_5_100

# For Make:
./compressor -s -o ../inputs/rMatGraph_J_5_100.binary ../inputs/rMatGraph_J_5_100
```

After an uncompressed graph has been converted to the binary format,
applications can be run on it by passing in the usual command-line flags, with
an additional `-b` flag. Note that the application will always load the binary
file using mmap.

```sh
# For Bazel:
$ bazel run //benchmarks/BFS/NonDeterministicBFS:BFS_main -- -s -b -src 10 ~/gbbs/inputs/rMatGraph_J_5_100.binary

# For Make:
$ ./BFS -s -b -src 10 ../../../inputs/rMatGraph_J_5_100.binary
```

Note that application performance will be affected if the file is not already
in the page-cache. We have found that using `-m` when the binary graph is backed
by SSD or disk results in a slow first-run, followed by fast subsequent runs.


Input Formats
-----------
We support the adjacency graph format used by the [Problem Based Benchmark
suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html)
and [Ligra](https://github.com/jshun/ligra).

The adjacency graph format starts with a sequence of offsets one for each
vertex, followed by a sequence of directed edges ordered by their source vertex.
The offset for a vertex i refers to the location of the start of a contiguous
block of out edges for vertex i in the sequence of edges. The block continues
until the offset of the next vertex, or the end if i is the last vertex. All
vertices and offsets are 0 based and represented in decimal. The specific format
is as follows:

```
AdjacencyGraph
<n>
<m>
<o0>
<o1>
...
<o(n-1)>
<e0>
<e1>
...
<e(m-1)>
```

This file is represented as plain text.

Weighted graphs are represented in the weighted adjacency graph format. The file
should start with the string "WeightedAdjacencyGraph". The m edge weights
should be stored after all of the edge targets in the .adj file.

**Using SNAP graphs**

Graphs from the [SNAP dataset
collection](https://snap.stanford.edu/data/index.html) are commonly used for
graph algorithm benchmarks. We provide a tool that converts the most common SNAP
graph format to the adjacency graph format that GBBS accepts. Usage example:
```sh
# Download a graph from the SNAP collection.
wget https://snap.stanford.edu/data/wiki-Vote.txt.gz
gzip --decompress ${PWD}/wiki-Vote.txt.gz
# Run the SNAP-to-adjacency-graph converter.
# Run with Bazel:
bazel run //utils:snap_converter -- -s -i ${PWD}/wiki-Vote.txt -o <output file>
# Or run with Make:
#   cd utils
#   make snap_converter
#   ./snap_converter -s -i <input file> -o <output file>
```

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

