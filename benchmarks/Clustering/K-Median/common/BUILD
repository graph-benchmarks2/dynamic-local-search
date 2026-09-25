licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "KMedianCommon",
    hdrs = ["KMedianCommon.h"],
    deps = [
        "//benchmarks/BFS/BFSdist:BFSdist",
        "//benchmarks/BFS/BFSdist:BFSdistSeq",
        "//benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaStepping",
        "//benchmarks/PositiveWeightSSSP/Dijkstra:Dijkstra",
        "//gbbs:gbbs",
        "@parlaylib//parlay:primitives",
        "@parlaylib//parlay:sequence",
    ],
)

cc_library(
    name = "StaticSwapTable",
    hdrs = ["StaticSwapTable.h"],
    deps = [
        ":KMedianCommon",
        "@parlaylib//parlay:sequence",
    ],
)

cc_library(
    name = "Initializers",
    hdrs = ["Initializers.h"],
    deps = [
        ":KMedianCommon",
        "//benchmarks/Clustering/K-Center/Gonzalez:Gonzalez",
        "//gbbs:gbbs",
        "@parlaylib//parlay:sequence",
    ],
)

cc_library(
    name = "RangeMaxTree",
    hdrs = ["RangeMaxTree.h"],
    deps = [
        "//gbbs:gbbs",
    ],
)

cc_binary(
    name = "RangeMaxTreeTest",
    srcs = ["RangeMaxTreeTest.cc"],
    deps = [
        ":RangeMaxTree",
    ],
)

cc_library(
    name = "DynamicSwapTable",
    hdrs = ["DynamicSwapTable.h"],
    deps = [
        ":RangeMaxTree",
        ":StaticSwapTable",
        "//gbbs:gbbs",
    ],
)

cc_binary(
    name = "DynamicSwapTableTest",
    srcs = [
        "DynamicSwapTableTest.cc",
        "DynamicSwapTableTestHelpers.h",
    ],
    deps = [
        ":DynamicSwapTable",
        ":KMedianCommon",
        ":StaticSwapTable",
        "//gbbs:gbbs",
    ],
)

cc_library(
    name = "DistanceDecrease",
    hdrs = ["DistanceDecrease.h"],
    deps = [
        "//gbbs:gbbs",
    ],
)

cc_library(
    name = "IncrementalAPSP",
    hdrs = ["IncrementalAPSP.h"],
    deps = [
        ":DistanceDecrease",
        ":KMedianCommon",
        "//benchmarks/BFS/BFSdist:BFSdistRepairSeq",
        "//benchmarks/BFS/BFSdist:BFSdistSeq",
        "//benchmarks/PositiveWeightSSSP/Dijkstra:Dijkstra",
        "//benchmarks/PositiveWeightSSSP/Dijkstra:DijkstraRepair",
        "//gbbs:gbbs",
    ],
)

cc_binary(
    name = "IncrementalAPSPTest",
    srcs = ["IncrementalAPSPTest.cc"],
    deps = [
        ":IncrementalAPSP",
        "//benchmarks/BFS/BFSdist:BFSdistSeq",
        "//benchmarks/PositiveWeightSSSP/Dijkstra:Dijkstra",
        "//gbbs:gbbs",
    ],
)
