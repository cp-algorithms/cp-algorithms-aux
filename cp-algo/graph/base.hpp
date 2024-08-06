#ifndef CP_ALGO_GRAPH_BASE_HPP
#define CP_ALGO_GRAPH_BASE_HPP
#include "edge_types.hpp"
#include "../data_structures/stack_union.hpp"
#include <ranges>
#include <vector>
namespace cp_algo::graph {
    using edge_index = int;
    enum type {directed = 0, undirected = 1};
    template<type _undirected, edge_type edge_t = edge_base>
    struct graph {
        static constexpr bool undirected = _undirected;
        graph(int n, int v0 = 0): v0(v0), adj(n) {}

        void add_edge(node_index u, edge_t e) {
            adj.push(u, size(edges));
            edges.push_back(e);
            if constexpr (undirected) {
                adj.push(e.to, size(edges));
            }
            edges.push_back(e.backedge(u));
        }
        void read_edges(node_index m) {
            adj.reserve(m);
            for(edge_index i = 0; i < m; i++) {
                auto [u, e] = edge_t::read(v0);
                add_edge(u, e);
            }
        }
        void call_adjacent(node_index v, auto &&callback, auto &&terminate) const {
            for(int sv = adj.head[v]; sv && !terminate(); sv = adj.next[sv]) {
                callback(adj.data[sv]);
            }
        }
        void call_adjacent(node_index v, auto &&callback) const {
            call_adjacent(v, callback, [](){return false;});
        }
        void call_edges(auto &&callback) const {
            for(edge_index e: edges_view()) {
                callback(e);
            }
        }
        auto nodes_view() const {
            return std::views::iota(0, n());
        }
        auto edges_view() const {
            return std::views::filter(
                std::views::iota(0, 2 * m()),
                [](edge_index e) {return !(e % 2);}
            );
        }
        auto const& incidence_lists() const {return adj;}
        edge_t const& edge(edge_index e) const {return edges[e];}
        node_index n() const {return adj.size();}
        edge_index m() const {return size(edges) / 2;}
    private:
        node_index v0;
        std::vector<edge_t> edges;
        data_structures::stack_union<edge_index> adj;
    };
}
#endif // CP_ALGO_GRAPH_BASE_HPP
