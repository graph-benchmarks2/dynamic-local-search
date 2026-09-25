licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "Dijkstra",
    hdrs = ["Dijkstra.h"],
    deps = [
        "//gbbs",
    ],
)

cc_library(
    name = "DijkstraRepair",
    hdrs = ["DijkstraRepair.h"],
    deps = [
        "//gbbs",
        "//benchmarks/PositiveWeightSSSP/common:PositiveWeightSSSPCommon",
    ],
)

cc_library(
    name = "DijkstraMultiSources",
    hdrs = ["DijkstraMultiSources.h"],
    deps = [
        "//gbbs",
    ],
)

cc_binary(
    name = "Dijkstra_main",
    srcs = ["Dijkstra.cc"],
    deps = [
        ":Dijkstra",
        "//gbbs",
    ],
)

cc_binary(
    name = "DijkstraRepair_main",
    srcs = ["DijkstraRepair.cc"],
    deps = [
        ":DijkstraRepair",
        "//gbbs",
    ],
)
