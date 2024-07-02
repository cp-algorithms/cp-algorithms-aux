#ifndef CP_ALGO_GRAPH_EULER_HPP
#define CP_ALGO_GRAPH_EULER_HPP
#include "base.hpp"
#include <algorithm>
#include <utility>
#include <vector>
namespace cp_algo::graph {
    template<typename graph>
    int euler_start(graph const& g) {
        std::vector<int> deg(g.n());
        constexpr bool undirected = graph::undirected;
        int res = 0;
        g.call_edges([&](int u, int e) {
            res = u;
            if constexpr (undirected) {
                deg[u] ^= 1;
            } else {
                deg[u]++;
                deg[g.edge(e).to]--;
            }
        });
        auto nodes = g.nodes_view();
        auto is_start = [&](int v) {return deg[v] > 0;};
        auto starts = std::ranges::count_if(nodes, is_start);
        if(starts > 1 + undirected) {
            res = -1;
        } else if(starts == 1 + undirected) {
            auto start = *std::ranges::find_if(nodes, is_start);
            res = deg[start] == 1 ? start : -1;
        }
        return res;
    }
    auto euler_trail(auto const& g) {
        int v0 = euler_start(g);
        std::vector<int> trail;
        if(~v0) {
            std::vector<bool> used(g.m());
            auto& adj = g.incidence_lists();
            auto head = adj.head;
            auto dfs = [&](auto &&self, int v) -> void {
                while(head[v]) {
                    int e = adj.data[std::exchange(head[v], adj.next[head[v]])];
                    if(!used[e / 2]) {
                        used[e / 2] = 1;
                        int u = g.edge(e).to;
                        self(self, u);
                        trail.push_back(e);
                    }
                }
            };
            dfs(dfs, v0);
            std::ranges::reverse(trail);
        }
        return std::pair{v0, trail};
    }
}
#endif // CP_ALGO_GRAPH_EULER_HPP
