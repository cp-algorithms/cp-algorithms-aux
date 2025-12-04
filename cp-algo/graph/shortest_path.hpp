#ifndef CP_ALGO_GRAPH_SHORTEST_PATH_HPP
#define CP_ALGO_GRAPH_SHORTEST_PATH_HPP
#include "base.hpp"
#include <algorithm>
#include <queue>
namespace cp_algo::graph {
    struct shortest_path_context {
        struct reverse {
            node_index u;
            edge_index e;
        };
        std::vector<uint64_t> dist;
        std::vector<reverse> pre;
        static constexpr uint64_t inf = 1e18;
        shortest_path_context(int n)
            : dist(n, inf), pre(n) {}
    };

    struct dijkstra_context : shortest_path_context {
        struct que_t {
            uint64_t dist;
            node_index v;
            bool operator<(que_t const& other) const {
                return dist > other.dist;
            }
        };
        std::priority_queue<que_t> pq;

        dijkstra_context(int n) : shortest_path_context(n) {}

        void push(node_index, node_index v) {
            pq.push({dist[v], v});
        }

        std::optional<node_index> next_node() {
            while (!empty(pq)) {
                auto [dv, v] = pq.top();
                pq.pop();
                if (dv == dist[v]) {
                    return v;
                }
            }
            return std::nullopt;
        }
    };

    struct spfa_context : shortest_path_context {
        std::queue<node_index> que;
        std::vector<char> flags;
        static constexpr char in_queue = 1;
        static constexpr char invalidated = 2;
        std::vector<std::basic_string<node_index>> dependents;

        spfa_context(int n) : shortest_path_context(n), flags(n), dependents(n) {}

        void push(node_index u, node_index v) {
            invalidate_subtree(v);
            dependents[u].push_back(v);
            flags[v] &= ~invalidated;
            if (!(flags[v] & in_queue)) {
                que.push(v);
                flags[v] |= in_queue;
            }
        }

        void invalidate_subtree(node_index v) {
            std::vector<node_index> to_invalidate = {v};
            while (!empty(to_invalidate)) {
                node_index u = to_invalidate.back();
                to_invalidate.pop_back();
                flags[u] |= invalidated;
                for (auto dep: dependents[u]) {
                    if (pre[dep].u == u) {
                        to_invalidate.push_back(dep);
                    }
                }
                dependents[u].clear();
            }
        }

        std::optional<node_index> next_node() {
            while (!que.empty()) {
                node_index v = que.front();
                que.pop();
                flags[v] &= ~in_queue;
                if (!(flags[v] & invalidated)) {
                    return v;
                }
            }
            return std::nullopt;
        }
    };

    template<typename Context, weighted_graph_type graph>
    Context sssp_impl(graph const& g, node_index s) {
        Context context(g.n());
        context.dist[s] = 0;
        context.push(s, s);
        while(auto ov = context.next_node()) {
            node_index v = *ov;
            for(auto e: g.outgoing(v)) {
                node_index u = g.edge(e).to;
                auto w = g.edge(e).w;
                if(context.dist[v] + w < context.dist[u]) {
                    context.dist[u] = context.dist[v] + w;
                    context.pre[u] = {v, e};
                    context.push(v, u);
                }
            }
        }
        return context;
    }

    template<weighted_graph_type graph>
    shortest_path_context single_source_shortest_path(graph const& g, node_index s) {
        bool negative_edges = false;
        for (auto e: g.edges()) {
            negative_edges |= e.w < 0;
        }
        if (negative_edges) {
            return sssp_impl<spfa_context>(g, s);
        } else {
            return sssp_impl<dijkstra_context>(g, s);
        }
    }

    template<weighted_graph_type graph>
    std::optional<std::pair<int64_t, std::vector<edge_index>>> shortest_path(graph const& g, node_index s, node_index t) {
        auto [dist, pre] = single_source_shortest_path(g, s);
        if (dist[t] == shortest_path_context::inf) {
            return std::nullopt;
        }
        std::vector<edge_index> path;
        node_index v = t;
        while(v != s) {
            path.push_back(pre[v].e);
            v = pre[v].u;
        }
        std::ranges::reverse(path);
        return {{dist[t], path}};
    }
}
#endif // CP_ALGO_GRAPH_SHORTEST_PATH_HPP
