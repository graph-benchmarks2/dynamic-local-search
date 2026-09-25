licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "RR",
    hdrs = ["RR.h"],
    deps = [
        "//gbbs",

        "//benchmarks/BFS/BFSdist:BFSdist",
        "//benchmarks/BFS/BFSdist:BFSdistRepair",
        "//benchmarks/BFS/BFSdist:BFSdistSeq",
        "//benchmarks/BFS/BFSdist:BFSdistRepairSeq",

        "//benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaStepping",
        "//benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaSteppingRepair",

        "//benchmarks/PositiveWeightSSSP/Dijkstra:Dijkstra",
        "//benchmarks/PositiveWeightSSSP/Dijkstra:DijkstraRepair",
    ],
)

cc_binary(
    name = "RR_main",
    srcs = ["RR.cc"],
    deps = [
        ":RR",
        "//benchmarks/Clustering/K-Center/common:DispatchMain",
    ],
)
