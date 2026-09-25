licenses(["notice"])

package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "DeltaStepping",
    hdrs = ["DeltaStepping.h"],
    deps = [
        "//gbbs",
        "//gbbs:bucket",
    ],
)

cc_library(
    name = "DeltaSteppingLight",
    hdrs = ["DeltaSteppingLight.h"],
    deps = [
        "//gbbs",
        "//gbbs:bucket",
    ],
)

cc_library(
    name = "DeltaSteppingRepair",
    hdrs = ["DeltaSteppingRepair.h"],
    deps = [
        "//benchmarks/PositiveWeightSSSP/common:PositiveWeightSSSPCommon",
        "//gbbs",
    ],
)

cc_library(
    name = "DeltaSteppingMultiSources",
    hdrs = ["DeltaSteppingMultiSources.h"],
    deps = [
        "//gbbs",
        "//gbbs:bucket",
    ],
)

cc_binary(
    name = "DeltaStepping_main",
    srcs = ["DeltaStepping.cc"],
    deps = [":DeltaStepping"],
)

cc_binary(
    name = "DeltaSteppingRepair_main",
    srcs = ["DeltaSteppingRepair.cc"],
    deps = [":DeltaSteppingRepair"],
)
