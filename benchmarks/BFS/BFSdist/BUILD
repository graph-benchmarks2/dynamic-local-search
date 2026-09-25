package(default_visibility = ["//visibility:public"])

cc_library(
    name = "BFSdist",
    hdrs = ["BFSdist.h"],
    deps = [
        "//gbbs",
    ],
)

cc_binary(
    name = "BFSdist_main",
    srcs = ["BFSdist.cc"],
    deps = [
        ":BFSdist",
        "//gbbs",
    ],
)

cc_library(
    name = "BFSdistRepair",
    hdrs = ["BFSdistRepair.h"],
    deps = [
        "//gbbs",
        "//benchmarks/PositiveWeightSSSP/common:PositiveWeightSSSPCommon",
    ],
)

cc_binary(
    name = "BFSdistRepair_main",
    srcs = ["BFSdistRepair.cc"],
    deps = [
        ":BFSdistRepair",
        "//gbbs",
    ],
)

cc_library(
    name = "BFSdistSeq",
    hdrs = ["BFSdistSeq.h"],
    deps = [
        "//gbbs",
    ],
)

cc_binary(
    name = "BFSdistSeq_main",
    srcs = ["BFSdistSeq.cc"],
    deps = [
        ":BFSdistSeq",
        "//gbbs",
    ],
)

cc_library(
    name = "BFSdistRepairSeq",
    hdrs = ["BFSdistRepairSeq.h"],
    deps = [
        "//gbbs",
        "//benchmarks/PositiveWeightSSSP/common:PositiveWeightSSSPCommon",
    ],
)

cc_binary(
    name = "BFSdistRepairSeq_main",
    srcs = ["BFSdistRepairSeq.cc"],
    deps = [
        ":BFSdistRepairSeq",
        "//gbbs",
    ],
)

cc_library(
    name = "BFSdistMultiSources",
    hdrs = ["BFSdistMultiSources.h"],
    deps = [
        "//gbbs",
    ],
)

cc_library(
    name = "BFSdistSeqMultiSources",
    hdrs = ["BFSdistSeqMultiSources.h"],
    deps = [
        "//gbbs",
    ],
)
