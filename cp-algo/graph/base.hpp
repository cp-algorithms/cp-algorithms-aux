#ifndef CP_ALGO_GRAPH_BASE_HPP
#define CP_ALGO_GRAPH_BASE_HPP
#include "../data_structures/stack_union.hpp"
#include <iostream>
#include <ranges>
#include <vector>
namespace cp_algo::graph {
    enum type {directed = 0, undirected = 1};
    template<type undirected>
    struct graph {
        graph(int n, int v0 = 0): n(n), v0(v0), adj(n) {}

        void add_edge(int u, int v) {
            adj.push(u - v0, size(to));
            to.push_back(v - v0);
            if constexpr (undirected) {
                adj.push(v - v0, size(to));
            }
            to.push_back(u - v0);
        }
        void read_edges(int m) {
            for(int i = 0; i < m; i++) {
                int u, v;
                std::cin >> u >> v;
                add_edge(u, v);
            }
        }
        auto nodes_view() const {
            return std::views::iota(0, n);
        }

        int n, v0;
        std::vector<int> to;
        data_structures::stack_union<int> adj;
    };
}
#endif // CP_ALGO_GRAPH_BASE_HPP
