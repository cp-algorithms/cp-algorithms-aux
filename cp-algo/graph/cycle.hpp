#ifndef CP_ALGO_GRAPH_CYCLE_HPP
#define CP_ALGO_GRAPH_CYCLE_HPP
#include "base.hpp"
#include <algorithm>
namespace cp_algo::graph {
    template<graph_type graph>
    std::vector<edge_index> find_cycle(graph const& g) {
        enum node_state { unvisited, visiting, visited };
        std::vector<node_state> state(g.n());
        std::vector<edge_index> cycle;
        bool closed = false;
        auto dfs = [&](this auto &&dfs, node_index v, edge_index ep = -1) -> bool {
            state[v] = visiting;
            for(auto e: g.outgoing(v)) {
                if constexpr (undirected_graph_type<graph>) {
                    if (ep == graph::opposite_idx(e)) {
                        continue;
                    }
                }
                node_index u = g.edge(e).to;
                if(state[u] == unvisited) {
                    if (dfs(u, e)) {
                        if (!closed) {
                            cycle.push_back(e);
                            closed |= g.edge(cycle[0]).to == v;
                        }
                        return true;
                    }
                } else if(state[u] == visiting) {
                    cycle = {e};
                    return true;
                }
            }
            state[v] = visited;
            return false;
        };
        for(node_index i: g.nodes()) {
            if(state[i] == unvisited && dfs(i)) {
                break;
            }
        }
        std::ranges::reverse(cycle);
        return cycle;
    }
}
#endif // CP_ALGO_GRAPH_CYCLE_HPP
