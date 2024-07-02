#ifndef CP_ALGO_GRAPH_CYCLE_HPP
#define CP_ALGO_GRAPH_CYCLE_HPP
#include "base.hpp"
#include <algorithm>
#include <vector>
namespace cp_algo::graph {
    std::vector<int> find_cycle(auto const& g) {
        std::vector<char> state(g.n());
        std::vector<int> cycle;
        auto dfs = [&](auto &&self, int v, int pe) -> bool {
            state[v] = 1;
            bool found = false;
            g.call_adjacent(v, [&](int e) {
                if(e / 2 != pe / 2) {
                    int u = g.edge(e).to;
                    if(state[u] == 0) {
                        if(self(self, u, e)) {
                            if(g.edge(cycle[0]).to != g.edge(cycle.back() ^ 1).to) {
                                cycle.push_back(e);
                            }
                            found = true;
                        }
                    } else if(state[u] == 1) {
                        cycle = {e};
                        found = true;
                    }
                }
            }, [&found](){return found;});
            state[v] = 2;
            return found;
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
