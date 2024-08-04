#ifndef CP_ALGO_GRAPH_SHORTEST_PATH_HPP
#define CP_ALGO_GRAPH_SHORTEST_PATH_HPP
#include "base.hpp"
#include <algorithm>
#include <queue>
namespace cp_algo::graph {
    template<type _, weighted_edge_type edge_t>
    auto dijkstra(graph<_, edge_t> const& g, int s) {
        static constexpr uint64_t inf = 1e18;
        std::vector<uint64_t> dist(g.n(), inf);
        std::vector<std::pair<node_index, edge_index>> pre(g.n());
        using que_t = std::pair<uint64_t, node_index>;
        std::priority_queue<que_t, std::vector<que_t>, std::greater<>> que;
        dist[s] = 0;
        que.push({0, s});
        while(!empty(que)) {
            auto [dv, v] = que.top();
            que.pop();
            if(dv != dist[v]) {
                continue;
            }
            g.call_adjacent(v, [&](auto e) {
                node_index u = g.edge(e).to;
                auto w = g.edge(e).w;
                if(dist[v] + w < dist[u]) {
                    pre[u] = {v, e};
                    dist[u] = dist[v] + w;
                    que.push({dist[u], u});
                }
            });
        }
        return std::pair{dist, pre};
    }

    template<type _, weighted_edge_type edge_t>
    auto shortest_path(graph<_, edge_t> const& g, int s, int t) {
        static constexpr uint64_t inf = 1e18;
        auto [dist, pre] = dijkstra(g, s);
        std::vector<std::pair<node_index, edge_index>> path;
        int64_t d = dist[t] == inf ? -1 : dist[t];
        while(d != -1 && t != s) {
            path.push_back(pre[t]);
            t = pre[t].first;
        }
        std::ranges::reverse(path);
        return std::pair{d, path};
    }
}
#endif // CP_ALGO_GRAPH_SHORTEST_PATH_HPP
