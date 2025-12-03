#ifndef CP_ALGO_GRAPH_2CC_HPP
#define CP_ALGO_GRAPH_2CC_HPP
#include "base.hpp"
#include "../structures/csr.hpp"
#include <algorithm>

namespace cp_algo::graph {
    template<undirected_graph_type graph>
    structures::csr<node_index> two_edge_connected_components(graph const& g) {
        structures::csr<node_index> components;
        components.reserve_data(g.n());
        enum node_state { unvisited, visited, blocked };
        std::vector<node_state> state(g.n());

        auto collect = [&](this auto &&collect, node_index v) -> void {
            state[v] = blocked;
            components.push(v);
            for (auto e : g.outgoing(v)) {
                node_index u = g.edge(e).to;
                if (state[u] != blocked) collect(u);
            }
        };

        big_vector<int> tin(g.n()), low(g.n());
        int timer = 0;
        auto dfs = [&](this auto &&dfs, node_index v, edge_index ep = -1) -> void {
            tin[v] = low[v] = timer++;
            state[v] = blocked;
            for (auto e: g.outgoing(v)) {
                if (ep == graph::opposite_idx(e)) {
                    continue;
                }
                node_index u = g.edge(e).to;
                if (state[u] == unvisited) {
                    dfs(u, e);
                    low[v] = std::min(low[v], low[u]);
                } else {
                    low[v] = std::min(low[v], tin[u]);
                }
            }
            if (low[v] == tin[v]) {
                components.new_row();
                collect(v);
            } else {
                state[v] = visited;
            }
        };

        for (auto v : g.nodes()) {
            if (state[v] == unvisited) {
                if (std::ranges::empty(g.outgoing(v))) {
                    components.new_row();
                    components.push(v);
                } else {
                    dfs(v);
                }
            }
        }
        return components;
    }
}

#endif // CP_ALGO_GRAPH_2CC_HPP