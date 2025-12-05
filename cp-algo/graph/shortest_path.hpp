#ifndef CP_ALGO_GRAPH_SHORTEST_PATH_HPP
#define CP_ALGO_GRAPH_SHORTEST_PATH_HPP
#include "base.hpp"
#include <algorithm>
#include <queue>
namespace cp_algo::graph {
    struct shortest_path_context {
        std::vector<int64_t> dist;
        std::vector<edge_index> pre;
        static constexpr int64_t inf = 1e18;
        shortest_path_context(int n)
            : dist(n, inf), pre(n) {}
    };

    struct dijkstra_context: shortest_path_context {
        struct que_t {
            int64_t dist;
            node_index v;
            bool operator<(que_t const& other) const {
                return dist > other.dist;
            }
        };
        std::priority_queue<que_t> pq;

        dijkstra_context(int n) : shortest_path_context(n) {}

        void push(node_index, edge_index, node_index v) {
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

    struct spfa_context: shortest_path_context {
        std::queue<node_index> que;
        std::vector<char> flags;
        static constexpr char in_queue = 1;
        static constexpr char invalidated = 2;

        spfa_context(int n) : shortest_path_context(n), flags(n) {}

        void push(node_index, edge_index, node_index v) {
            if (!(flags[v] & in_queue)) {
                que.push(v);
                flags[v] |= in_queue;
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

    struct deep_spfa_context: spfa_context {
        struct traverse_edge {
            edge_index e;
            node_index v;
        };
        std::vector<std::basic_string<traverse_edge>> dependents;

        deep_spfa_context(int n) : spfa_context(n), dependents(n) {}

        void push(node_index u, edge_index e, node_index v) {
            invalidate_subtree(v);
            dependents[u].push_back({e, v});
            flags[v] &= ~invalidated;
            spfa_context::push(u, e, v);
        }

        void invalidate_subtree(node_index v) {
            std::vector<node_index> to_invalidate = {v};
            while (!empty(to_invalidate)) {
                node_index u = to_invalidate.back();
                to_invalidate.pop_back();
                flags[u] |= invalidated;
                flags[u] &= ~in_queue;
                for (auto [e, v]: dependents[u]) {
                    if (pre[v] == e) {
                        to_invalidate.push_back(v);
                    }
                }
                dependents[u].clear();
            }
        }
    };

    template<typename Context, weighted_graph_type graph>
    Context sssp_impl(graph const& g, node_index s) {
        Context context(g.n());
        context.dist[s] = 0;
        context.pre[s] = -1;
        context.push(s, -1, s);
        while(auto ov = context.next_node()) {
            node_index v = *ov;
            for(auto e: g.outgoing(v)) {
                node_index u = g.edge(e).traverse(v);
                auto w = g.edge(e).w;
                if(context.dist[v] + w < context.dist[u]) {
                    context.dist[u] = context.dist[v] + w;
                    context.pre[u] = e;
                    context.push(v, e, u);
                }
            }
        }
        return context;
    }

    template<weighted_graph_type graph>
    shortest_path_context dijkstra(graph const& g, node_index s) {
        return sssp_impl<dijkstra_context>(g, s);
    }
    template<weighted_graph_type graph>
    shortest_path_context spfa(graph const& g, node_index s) {
        return sssp_impl<spfa_context>(g, s);
    }
    template<weighted_graph_type graph>
    shortest_path_context deep_spfa(graph const& g, node_index s) {
        return sssp_impl<deep_spfa_context>(g, s);
    }
    
    template<weighted_graph_type graph>
    shortest_path_context single_source_shortest_path(graph const& g, node_index s) {
        bool negative_edges = false;
        for (auto e: g.edges()) {
            negative_edges |= e.w < 0;
        }
        return negative_edges ? deep_spfa(g, s) : dijkstra(g, s);
    }

    std::vector<edge_index> recover_path(auto const& g, auto const& pre, node_index s, node_index t) {
        std::vector<edge_index> path;
        node_index v = t;
        while(v != s) {
            path.push_back(pre[v]);
            v = g.edge(pre[v]).traverse(v);
        }
        std::ranges::reverse(path);
        return path;
    }

    template<weighted_graph_type graph>
    std::optional<std::pair<int64_t, std::vector<edge_index>>> shortest_path(graph const& g, node_index s, node_index t) {
        auto [dist, pre] = single_source_shortest_path(g, s);
        if (dist[t] == shortest_path_context::inf) {
            return std::nullopt;
        }
        return {{dist[t], recover_path(g, pre, s, t)}};
    }
}
#endif // CP_ALGO_GRAPH_SHORTEST_PATH_HPP
