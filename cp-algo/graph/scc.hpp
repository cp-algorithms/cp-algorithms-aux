#ifndef CP_ALGO_GRAPH_SCC_HPP
#define CP_ALGO_GRAPH_SCC_HPP
#include "base.hpp"
#include "../structures/csr.hpp"
#include <algorithm>
namespace cp_algo::graph {
    // Tarjan's algorithm for Strongly Connected Components
    // returns components in reverse topological order
    template<digraph_type graph>
    structures::csr<node_index> scc(graph const& g) {
        structures::csr<node_index> components;
        components.reserve_data(g.n());
        enum node_state { unvisited, visited, collected };
        std::vector<node_state> state(g.n());
        auto collect = [&](this auto &&collect, node_index v) -> void {
            state[v] = collected;
            components.push(v);
            for(auto e: g.outgoing(v)) {
                node_index u = g.edge(e).to;
                if (state[u] != collected) {
                    collect(u);
                }
            }
        };
        big_vector<int> tin(g.n()), low(g.n());
        int timer = 0;
        auto dfs = [&](this auto &&dfs, node_index v) -> void {
            state[v] = visited;
            tin[v] = low[v] = timer++;
            for(auto e: g.outgoing(v)) {
                node_index u = g.edge(e).to;
                if (state[u] == unvisited) {
                    dfs(u);
                    low[v] = std::min(low[v], low[u]);
                } else if (state[u] == visited) {
                    low[v] = std::min(low[v], tin[u]);
                }
            }
            if (low[v] == tin[v]) {
                components.new_row();
                collect(v);
            }
        };
        for (auto v: g.nodes()) {
            if (state[v] == unvisited) {
                dfs(v);
            }
        }
        return components;
    }
}
#endif // CP_ALGO_GRAPH_SCC_HPP
