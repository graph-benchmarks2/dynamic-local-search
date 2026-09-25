licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "Gonzalez",
    hdrs = ["Gonzalez.h"],
    deps = [
        "//benchmarks/BFS/BFSdist:BFSdist",
        "//benchmarks/BFS/BFSdist:BFSdistRepair",
        "//benchmarks/BFS/BFSdist:BFSdistSeq",
        "//benchmarks/BFS/BFSdist:BFSdistRepairSeq",

        "//benchmarks/Clustering/K-Center/Gonzalez/common:KCenterResult",

        "//benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaStepping",
        "//benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaSteppingRepair",

        "//benchmarks/PositiveWeightSSSP/Dijkstra:Dijkstra",
        "//benchmarks/PositiveWeightSSSP/Dijkstra:DijkstraRepair",

        "//gbbs:bucket",
        "//gbbs/helpers:dyn_arr",
        "//gbbs:graph",

        "@parlaylib//parlay:parallel",
        "@parlaylib//parlay:primitives",
        "@parlaylib//parlay:random",
        "@parlaylib//parlay:sequence",
        "@parlaylib//parlay:utilities",
    ],
)

cc_binary(
    name = "Gonzalez_main",
    srcs = ["Gonzalez.cc"],
    deps = [
        ":Gonzalez",
        "//benchmarks/Clustering/K-Center/common:DispatchMain",
        "//benchmarks/Clustering/K-Center/common:kcenter_quality",
    ],
)
