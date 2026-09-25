#pragma once

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

namespace gbbs {
namespace kcenter_common {

enum class GraphWeightFormat {
  kUnweighted,
  kWeighted
};

inline const char* graph_weight_format_name(GraphWeightFormat fmt) {
  return (fmt == GraphWeightFormat::kWeighted) ? "weighted" : "unweighted";
}

// Detect weighted/unweighted by inspecting the actual graph file header.
// Currently supports only plain text adjacency-graph files.
//
// Supported headers:
//   AdjacencyGraph
//   WeightedAdjacencyGraph
//
// For compressed (-c) and binary (-b) graphs, we currently do not have a clean
// pre-load format check here, so we fail explicitly rather than guessing.
inline GraphWeightFormat detect_graph_weight_format_or_die(
    const char* filename, bool compressed, bool binary) {
  if (filename == nullptr) {
    std::cout << "ERROR: No input graph file was provided." << std::endl;
    std::exit(-1);
  }

  if (compressed) {
    std::cout
        << "ERROR: Automatic weighted/unweighted format detection is not yet "
           "implemented for compressed graphs (-c).\n"
        << "Please use a plain text graph file for now."
        << std::endl;
    std::exit(-1);
  }

  if (binary) {
    std::cout
        << "ERROR: Automatic weighted/unweighted format detection is not yet "
           "implemented for binary graphs (-b).\n"
        << "Please use a plain text graph file for now."
        << std::endl;
    std::exit(-1);
  }

  std::ifstream in(filename);
  if (!in.is_open()) {
    std::cout << "ERROR: Unable to open input graph file: " << filename
              << std::endl;
    std::exit(-1);
  }

  std::string header;
  in >> header;

  if (header.empty()) {
    std::cout << "ERROR: Input graph file appears to be empty: " << filename
              << std::endl;
    std::exit(-1);
  }

  if (header == "AdjacencyGraph") {
    return GraphWeightFormat::kUnweighted;
  }

  if (header == "WeightedAdjacencyGraph") {
    return GraphWeightFormat::kWeighted;
  }

  std::cout
      << "ERROR: Unrecognized graph header in input file.\n"
      << "Input file: " << filename << "\n"
      << "Header token: " << header << "\n"
      << "Expected one of: AdjacencyGraph, WeightedAdjacencyGraph"
      << std::endl;
  std::exit(-1);
}

}  // namespace kcenter_common
}  // namespace gbbs
