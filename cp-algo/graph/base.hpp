#ifndef CP_ALGO_GRAPH_BASE_HPP
#define CP_ALGO_GRAPH_BASE_HPP
#include "../data_structures/stack_union.hpp"
#include <iostream>
#include <ranges>
#include <vector>
namespace cp_algo::graph {
    using node_index = int;
    using edge_index = int;
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

    enum type {directed = 0, undirected = 1};
    template<type _undirected, edge_type edge_t = edge_base>
    struct graph {
        static constexpr bool undirected = _undirected;
        graph(int n, int v0 = 0): v0(v0), _adj(n) {}

        void add_edge(node_index u, edge_t e) {
            _adj.push(u, size(edges));
            edges.push_back(e);
            if constexpr (undirected) {
                _adj.push(e.to, size(edges));
            }
            edges.push_back(e.backedge(u));
        }
        void read_edges(node_index m) {
            for(edge_index i = 0; i < m; i++) {
                auto [u, e] = edge_t::read(v0);
                add_edge(u, e);
            }
        }
        void call_adjacent(node_index v, auto &&callback, auto &&terminate) const {
            for(int sv = _adj.head[v]; sv && !terminate(); sv = _adj.next[sv]) {
                callback(_adj.data[sv]);
            }
        }
        void call_adjacent(node_index v, auto &&callback) const {
            call_adjacent(v, callback, [](){return false;});
        }
        void call_edges(auto &&callback) const {
            for(node_index v: nodes_view()) {
                call_adjacent(v, [&](edge_index e) {callback(v, e);});
            }
        }
        auto nodes_view() const {
            return std::views::iota(0, n());
        }
        auto const& incidence_lists() const {return _adj;}
        edge_t const& edge(edge_index e) const {return edges[e];}
        size_t n() const {return _adj.size();}
        size_t m() const {return size(edges) / 2;}
    private:
        node_index v0;
        std::vector<edge_t> edges;
        data_structures::stack_union<edge_index> _adj;
    };
}
#endif // CP_ALGO_GRAPH_BASE_HPP
