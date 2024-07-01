#ifndef CP_ALGO_GRAPH_CYCLE_HPP
#define CP_ALGO_GRAPH_CYCLE_HPP
#include "base.hpp"
#include <algorithm>
#include <vector>
namespace cp_algo::graph {
    std::vector<size_t> find_cycle(auto const& g) {
        std::vector<char> state(g.n);
        std::vector<size_t> cycle;
        auto dfs = [&](auto self, size_t v, size_t pe) -> bool {
            state[v] = 1;
            auto gen = g.adjacent_generator(v);
            for(size_t e = gen(); ~e; e = gen()) {
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
        for(size_t i: g.nodes_view()) {
            if(!state[i] && dfs(dfs, i, -1)) {
                break;
            }
        }
        std::ranges::reverse(cycle);
        return cycle;
    }
}
#endif // CP_ALGO_GRAPH_CYCLE_HPP
