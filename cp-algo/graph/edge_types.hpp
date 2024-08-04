#ifndef CP_ALGO_GRAPH_EDGE_TYPES_HPP
#define CP_ALGO_GRAPH_EDGE_TYPES_HPP
#include <iostream>
#include <cstdint>
namespace cp_algo::graph {
    using node_index = int;
    struct edge_base {
        node_index to;

        edge_base() {}
        edge_base(node_index v): to(v) {}

        static auto read(node_index v0 = 0) {
            node_index u, v;
            std::cin >> u >> v;
            return std::pair{u - v0, edge_base(v - v0)};
        }

        edge_base backedge(int from) const {
            return {from};
        }
    };

    struct weighted_edge: edge_base {
        int64_t w;

        weighted_edge() {}
        weighted_edge(node_index v, int64_t w): edge_base(v), w(w) {}

        static auto read(node_index v0 = 0) {
            node_index u, v;
            int64_t w;
            std::cin >> u >> v >> w;
            return std::pair{u - v0, weighted_edge{v - v0, w}};
        }

        weighted_edge backedge(node_index from) const {
            return {from, w};
        }
    };

    template<typename edge>
    concept edge_type = std::is_base_of_v<edge_base, edge>;
    template<typename edge>
    concept weighted_edge_type = std::is_base_of_v<weighted_edge, edge>;
}
#endif // CP_ALGO_GRAPH_EDGE_TYPES_HPP
