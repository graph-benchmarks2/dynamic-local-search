#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "gbbs/gbbs.h"
#include "parlay/sequence.h"

namespace gbbs {
namespace kmedian {

template <class W>
struct DynamicEdgeInsertion {
    uintE u;
    uintE v;
    W w;
};

inline uint64_t DynamicUndirectedEdgeKey(uintE u, uintE v) {
    if (u > v) {
        std::swap(u, v);
    }

    return
        (static_cast<uint64_t>(u) << 32)
        | static_cast<uint64_t>(v);
}

template <class Graph>
parlay::sequence<typename Graph::edge>
ExtractDynamicUniqueUndirectedEdges(Graph& G) {
    using W = typename Graph::weight_type;
    using edge = typename Graph::edge;

    std::vector<edge> edges;
    edges.reserve(G.m / 2 + 1);

    for (uintE u = 0; u < G.n; ++u) {
        auto map_f =
            [&](const uintE& src, const uintE& v, const W& w) {

                if (src < v) {
                    edges.emplace_back(src, v, w);
                }
            };

        G.get_vertex(u).out_neighbors().map(
            map_f,
            false);
    }

    return parlay::sequence<edge>(
        edges.begin(),
        edges.end());
}

template <class Graph>
parlay::sequence<DynamicEdgeInsertion<typename Graph::weight_type>>
ReadDynamicInsertions(
    const std::string& filename,
    size_t n) {

    using W = typename Graph::weight_type;
    using insertion = DynamicEdgeInsertion<W>;

    std::ifstream in(filename);

    if (!in) {
        std::cout << "ERROR: Could not open insertion file: "
                  << filename
                  << std::endl;
        std::exit(-1);
    }

    std::vector<insertion> updates;

    while (true) {
        uint64_t u_raw;
        uint64_t v_raw;

        if (!(in >> u_raw >> v_raw)) {
            break;
        }

        if (u_raw >= n || v_raw >= n) {
            std::cout << "ERROR: Insertion endpoint outside graph: "
                      << u_raw << " " << v_raw
                      << std::endl;
            std::exit(-1);
        }

        if (u_raw == v_raw) {
            std::cout << "ERROR: Self-loop insertion is not supported: "
                      << u_raw << " " << v_raw
                      << std::endl;
            std::exit(-1);
        }

        const uintE u =
            static_cast<uintE>(u_raw);

        const uintE v =
            static_cast<uintE>(v_raw);

        if constexpr (std::is_same<W, gbbs::empty>::value) {
            updates.push_back(
                insertion{u, v, W{}});
        } else {
            double w_raw;

            if (!(in >> w_raw)) {
                std::cout << "ERROR: Weighted insertion file must use 'u v w'."
                          << std::endl;
                std::exit(-1);
            }

            if (w_raw <= 0.0) {
                std::cout << "ERROR: Edge weights must be positive."
                          << std::endl;
                std::exit(-1);
            }

            updates.push_back(
                insertion{
                    u,
                    v,
                    static_cast<W>(w_raw)});
        }
    }

    return parlay::sequence<insertion>(
        updates.begin(),
        updates.end());
}

template <class W>
using DynamicTempGraph =
    gbbs::symmetric_graph<
        gbbs::symmetric_vertex,
        W>;

template <class W>
DynamicTempGraph<W> BuildDynamicGraphFromEdges(
    const parlay::sequence<std::tuple<uintE, uintE, W>>& edges,
    size_t n) {

    return DynamicTempGraph<W>::from_edges(
        edges,
        n);
}

}  // namespace kmedian
}  // namespace gbbs
