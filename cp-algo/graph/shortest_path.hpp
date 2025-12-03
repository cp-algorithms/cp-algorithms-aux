#ifndef CP_ALGO_GRAPH_SHORTEST_PATH_HPP
#define CP_ALGO_GRAPH_SHORTEST_PATH_HPP
#include "base.hpp"
#include <algorithm>
#include <queue>
namespace cp_algo::graph {
    template<weighted_graph_type graph>
    auto dijkstra(graph const& g, int s) {
        static constexpr uint64_t inf = 1e18;
        std::vector<uint64_t> dist(g.n(), inf);
        struct reverse {
            node_index u;
            edge_index e;
        };
        std::vector<reverse> pre(g.n());
        using que_t = std::pair<uint64_t, node_index>;
        std::priority_queue<que_t, std::vector<que_t>, std::greater<>> que;
        dist[s] = 0;
        que.push({0, s});
        std::vector<bool> visited(g.n());
        while(!empty(que)) {
            auto [dv, v] = que.top();
            que.pop();
            if(dv != dist[v]) {
                continue;
            }
            for(auto e: g.outgoing(v)) {
                node_index u = g.edge(e).to;
                auto w = g.edge(e).w;
                if(dist[v] + w < dist[u]) {
                    pre[u] = {v, e};
                    dist[u] = dist[v] + w;
                    que.push({dist[u], u});
                }
            }
        }
        return std::pair{dist, pre};
    }

    template<weighted_graph_type graph>
    std::optional<std::pair<int64_t, std::vector<edge_index>>> shortest_path(graph const& g, int s, int t) {
        static constexpr uint64_t inf = 1e18;
        auto [dist, pre] = dijkstra(g, s);
        if (dist[t] == inf) {
            return std::nullopt;
        }
        std::vector<edge_index> path;
        int64_t d = dist[t];
        while(t != s) {
            path.push_back(pre[t].e);
            t = pre[t].u;
        }
        std::ranges::reverse(path);
        return {{d, path}};
    }
}
#endif // CP_ALGO_GRAPH_SHORTEST_PATH_HPP
