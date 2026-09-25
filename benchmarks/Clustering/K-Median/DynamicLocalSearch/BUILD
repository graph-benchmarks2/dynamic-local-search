licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "DynamicLocalSearch",
    hdrs = ["DynamicLocalSearch.h"],
    deps = [
        "//benchmarks/Clustering/K-Median/common:DistanceDecrease",
        "//benchmarks/Clustering/K-Median/common:DynamicSwapTable",
        "//benchmarks/Clustering/K-Median/common:KMedianCommon",
        "//gbbs:gbbs",
    ],
)

cc_binary(
    name = "DynamicLocalSearchTest",
    srcs = ["DynamicLocalSearchTest.cc"],
    deps = [
        ":DynamicLocalSearch",
        "//benchmarks/Clustering/K-Median/common:DynamicSwapTable",
        "//benchmarks/Clustering/K-Median/common:KMedianCommon",
        "//gbbs:gbbs",
    ],
)

cc_library(
    name = "DynamicGraphIO",
    hdrs = ["DynamicGraphIO.h"],
    deps = ["//gbbs:gbbs"],
)

cc_library(
    name = "DynamicGraphDriverTypes",
    hdrs = ["DynamicGraphDriverTypes.h"],
    deps = [
        "//benchmarks/Clustering/K-Median/common:KMedianCommon",
        "//gbbs:gbbs",
    ],
)

cc_library(
    name = "DynamicGraphReporting",
    hdrs = ["DynamicGraphReporting.h"],
    deps = [":DynamicGraphDriverTypes"],
)

cc_library(
    name = "DynamicGraphDriver",
    hdrs = ["DynamicGraphDriver.h"],
    deps = [
        ":DynamicGraphDriverTypes",
        ":DynamicGraphIO",
        ":DynamicGraphReporting",
        ":DynamicLocalSearch",
        "//benchmarks/Clustering/K-Median/StaticLocalSearch:StaticLocalSearch",
        "//benchmarks/Clustering/K-Median/common:DynamicSwapTable",
        "//benchmarks/Clustering/K-Median/common:IncrementalAPSP",
        "//benchmarks/Clustering/K-Median/common:KMedianCommon",
        "//benchmarks/Clustering/K-Median/common:StaticSwapTable",
        "//gbbs:gbbs",
    ],
)

cc_binary(
    name = "DynamicLocalSearch_main",
    srcs = ["DynamicLocalSearchMain.cc"],
    deps = [
        ":DynamicGraphDriver",
        "//benchmarks/Clustering/K-Center/common:DispatchMain",
        "//benchmarks/Clustering/K-Median/common:Initializers",
        "//gbbs:gbbs",
    ],
)
