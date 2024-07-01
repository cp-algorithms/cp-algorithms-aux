#ifndef CP_ALGO_GRAPH_BASE_HPP
#define CP_ALGO_GRAPH_BASE_HPP
#include <iostream>
#include <ranges>
#include <vector>
namespace cp_algo::graph {
    enum type {directed = 0, undirected = 1};
    template<type undirected>
    struct graph {
        graph(size_t n, size_t v0 = 0): n(n), v0(v0), adj(n) {}

        void add_edge(size_t u, size_t v) {
            adj[u - v0].push_back(size(to));
            to.push_back(v - v0);
            if constexpr (undirected) {
                adj[v - v0].push_back(size(to));
            }
            to.push_back(u - v0);
        }
        void read_edges(size_t m) {
            for(size_t i = 0; i < m; i++) {
                int u, v;
                std::cin >> u >> v;
                add_edge(u, v);
            }
        }
        auto adjacent_generator(size_t v) const {
            return [&adj = adj[v], idx = 0]() mutable {
                return idx < ssize(adj) ? adj[idx++] : -1;
            };
        }
        auto nodes_view() const {
            return std::views::iota((size_t)0, n);
        }

        size_t n;
        size_t v0;
        std::vector<size_t> to;
        std::vector<std::vector<size_t>> adj;
    };
}
#endif // CP_ALGO_GRAPH_BASE_HPP
