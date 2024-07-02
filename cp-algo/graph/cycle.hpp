#ifndef CP_ALGO_GRAPH_CYCLE_HPP
#define CP_ALGO_GRAPH_CYCLE_HPP
#include "base.hpp"
#include <algorithm>
#include <vector>
namespace cp_algo::graph {
    std::vector<int> find_cycle(auto const& g) {
        std::vector<char> state(g.n);
        std::vector<int> cycle;
        auto dfs = [&](auto &&self, int v, int pe) -> bool {
            state[v] = 1;
            for(int sv = g.adj.head[v]; sv; sv = g.adj.next[sv]) {
                int e = g.adj.data[sv];
                if(e / 2 != pe / 2) {
                    auto u = g.to[e];
                    if(state[u] == 0) {
                        if(self(self, u, e)) {
                            if(g.to[cycle[0]] != g.to[cycle.back() ^ 1]) {
                                cycle.push_back(e);
                            }
                            return true;
                        }
                    } else if(state[u] == 1) {
                        cycle = {e};
                        return true;
                    }
                }
            }
            state[v] = 2;
            return false;
        };
        for(int i: g.nodes_view()) {
            if(!state[i] && dfs(dfs, i, -2)) {
                break;
            }
        }
        std::ranges::reverse(cycle);
        return cycle;
    }
}
#endif // CP_ALGO_GRAPH_CYCLE_HPP
