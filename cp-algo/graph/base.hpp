#ifndef CP_ALGO_GRAPH_BASE_HPP
#define CP_ALGO_GRAPH_BASE_HPP
#include "../data_structures/stack_union.hpp"
#include <iostream>
#include <ranges>
#include <vector>
namespace cp_algo::graph {
    struct edge_base {
        int to;

        edge_base() {}
        edge_base(int v): to(v) {}

        static auto read(int v0 = 0) {
            int u, v;
            std::cin >> u >> v;
            return std::pair{u - v0, edge_base(v - v0)};
        }

        edge_base backedge(int from) const {
            return {from};
        }
    };
    template<typename edge_info>
    concept edge_type = std::is_base_of_v<edge_base, edge_info>;

    enum type {directed = 0, undirected = 1};
    template<type _undirected, edge_type edge_info = edge_base>
    struct graph {
        static constexpr bool undirected = _undirected;
        graph(int n, int v0 = 0): _n(n), _m(0), v0(v0), _adj(n) {}

        void add_edge(int u, edge_info e) {
            _m++;
            _adj.push(u, size(edges));
            edges.push_back(e);
            if constexpr (undirected) {
                _adj.push(e.to, size(edges));
            }
            edges.push_back(e.backedge(u));
        }
        void read_edges(int m) {
            for(int i = 0; i < m; i++) {
                auto [u, e] = edge_info::read(v0);
                add_edge(u, e);
            }
        }
        void call_adjacent(int v, auto &&callback, auto &&terminate) const {
            for(int sv = _adj.head[v]; sv && !terminate(); sv = _adj.next[sv]) {
                callback(_adj.data[sv]);
            }
        }
        void call_adjacent(int v, auto &&callback) const {
            call_adjacent(v, callback, [](){return false;});
        }
        void call_edges(auto &&callback) const {
            for(int v: nodes_view()) {
                call_adjacent(v, [&](int e) {callback(v, e);});
            }
        }
        auto nodes_view() const {
            return std::views::iota(0, _n);
        }
        auto const& incidence_lists() const {return _adj;}
        auto const& edge(int e) const {return edges[e];}
        int n() const {return _n;}
        int m() const {return _m;}
    private:
        int _n, _m, v0;
        std::vector<edge_info> edges;
        data_structures::stack_union<int> _adj;
    };
}
#endif // CP_ALGO_GRAPH_BASE_HPP
