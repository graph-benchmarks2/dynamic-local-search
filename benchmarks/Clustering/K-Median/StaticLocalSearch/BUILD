licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "StaticLocalSearch",
    hdrs = ["StaticLocalSearch.h"],
    deps = [
        "//benchmarks/Clustering/K-Median/common:KMedianCommon",
        "//benchmarks/Clustering/K-Median/common:StaticSwapTable",
        "//gbbs:gbbs",
    ],
)

cc_binary(
    name = "StaticLocalSearch_main",
    srcs = ["StaticLocalSearch.cc"],
    deps = [
        ":StaticLocalSearch",
        "//benchmarks/Clustering/K-Center/common:DispatchMain",
        "//benchmarks/Clustering/K-Median/common:IncrementalAPSP",
        "//benchmarks/Clustering/K-Median/common:Initializers",
    ],
)

cc_binary(
    name = "WarmStartStaticLocalSearch_main",
    srcs = [
        "WarmStartStaticLocalSearch.cc",
        "WarmStartDriver.h",
    ],
    deps = [
        ":StaticLocalSearch",
        "//benchmarks/Clustering/K-Median/common:IncrementalAPSP",
        "//benchmarks/Clustering/K-Median/common:Initializers",
        "//benchmarks/Clustering/K-Center/common:DispatchMain",
        "//gbbs:gbbs",
    ],
)
