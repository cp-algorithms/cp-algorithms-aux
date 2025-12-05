#ifndef CP_ALGO_GRAPH_EDGE_TYPES_HPP
#define CP_ALGO_GRAPH_EDGE_TYPES_HPP
#include <iostream>
#include <cstdint>
namespace cp_algo::graph {
    using node_index = int;
    struct edge_base {
        int xor_nodes;

        edge_base() {}
        edge_base(node_index from, node_index to): xor_nodes(from ^ to) {}

        // Given one endpoint, return the other
        node_index traverse(node_index from) const {
            return xor_nodes ^ from;
        }

        static auto read(node_index v0 = 0) {
            node_index u, v;
            std::cin >> u >> v;
            u -= v0;
            v -= v0;
            return std::pair{u, edge_base(u, v)};
        }
    };

    struct weighted_edge: edge_base {
        int64_t w;

        weighted_edge() {}
        weighted_edge(node_index from, node_index to, int64_t w): edge_base(from, to), w(w) {}

        static auto read(node_index v0 = 0) {
            auto [u, e] = edge_base::read(v0);
            int64_t w;
            std::cin >> w;
            return std::pair{u, weighted_edge(u, e.traverse(u), w)};
        }
    };

    template<typename edge>
    concept edge_type = std::is_base_of_v<edge_base, edge>;
    template<typename edge>
    concept weighted_edge_type = std::is_base_of_v<weighted_edge, edge>;
}
#endif // CP_ALGO_GRAPH_EDGE_TYPES_HPP
