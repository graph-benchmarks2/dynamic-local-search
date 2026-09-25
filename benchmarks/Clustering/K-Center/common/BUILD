licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "GraphFormatGuard",
    hdrs = ["GraphFormatGuard.h"],
)

cc_library(
    name = "DispatchMain",
    hdrs = ["DispatchMain.h"],
    deps = [
        ":GraphFormatGuard",
        "//gbbs",
    ],
)

cc_library(
    name = "kcenter_quality",
    hdrs = ["KCenterQuality.h"],
    deps = [
        "//benchmarks/BFS/BFSdist:BFSdist",
        "//benchmarks/BFS/BFSdist:BFSdistSeq",
        "//benchmarks/PositiveWeightSSSP/DeltaStepping:DeltaStepping",
        "//benchmarks/PositiveWeightSSSP/Dijkstra:Dijkstra",
        "//gbbs:gbbs",
    ],
)
