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
    template<type undirected, edge_type edge_info = edge_base>
    struct graph {
        graph(int n, int v0 = 0): _n(n), m(0), v0(v0), adj(n) {}

        void add_edge(int u, edge_info e) {
            m++;
            adj.push(u, size(edges));
            edges.push_back(e);
            if constexpr (undirected) {
                adj.push(e.to, size(edges));
            }
            edges.push_back(e.backedge(u));
        }
        void read_edges(int m) {
            for(int i = 0; i < m; i++) {
                auto [u, e] = edge_info::read(v0);
                add_edge(u, e);
            }
        }
        void call_adjacent(int v, auto &&callback, auto &&terminate = [](){return false;}) const {
            for(int sv = adj.head[v]; sv && !terminate(); sv = adj.next[sv]) {
                callback(adj.data[sv]);
            }
        }
        int deg(int v) const {
            int ans = 0;
            call_adjacent([&ans](){ans++;});
            return ans;
        }
        auto nodes_view() const {
            return std::views::iota(0, _n);
        }
        edge_info& edge(int e) {return edges[e];}
        edge_info const& edge(int e) const {return edges[e];}
        int n() const {return _n;}
    private:
        int _n, m, v0;
        std::vector<edge_info> edges;
        data_structures::stack_union<int> adj;
    };
}
#endif // CP_ALGO_GRAPH_BASE_HPP
