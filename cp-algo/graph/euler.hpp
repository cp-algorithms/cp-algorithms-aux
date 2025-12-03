#ifndef CP_ALGO_GRAPH_EULER_HPP
#define CP_ALGO_GRAPH_EULER_HPP
#include "base.hpp"
#include <algorithm>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>
#include <cassert>
namespace cp_algo::graph {
    template<graph_type graph>
    std::optional<node_index> euler_start(graph const& g) {
        std::vector<int> deg(g.n());
        std::optional<node_index> default_start = 0;
        for(auto v: g.nodes()) {
            for(auto e: g.outgoing(v)) {
                deg[v]++;
                default_start = v;
                if constexpr (digraph_type<graph>) {
                    deg[g.edge(e).to]--;
                }
            }
        }
        if constexpr (undirected_graph_type<graph>) {
            for(auto &it: deg) {
                it %= 2;
            }
        }
        auto is_start = [&](int v) {return deg[v] > 0;};
        auto starts = std::ranges::count_if(g.nodes(), is_start);
        auto need_starts = undirected_graph_type<graph> ? 2 : 1;
        if(starts > need_starts) {
            return std::nullopt;
        } else if(starts == need_starts) {
            auto start = *std::ranges::find_if(g.nodes(), is_start);
            if (deg[start] == 1) {
                return start;
            } else {
                return std::nullopt;
            }
        }
        return default_start;
    }
    // Try finding a trail starting from v0
    // may be partial if graph is not Eulerian or disconnected
    template<graph_type graph>
    std::vector<edge_index> try_euler_trail(graph const& g, node_index v0) {
        std::vector<edge_index> trail;
        std::vector<bool> used(g.m());
        auto head = g.nodes() | std::views::transform([&](auto v) {
            return begin(g.outgoing(v));
        }) | std::ranges::to<std::vector>();
        auto dfs = [&](this auto &&dfs, int v) -> void {
            while (head[v] != end(g.outgoing(v))) {
                auto e = *head[v]++;
                if(!used[graph::canonical_idx(e)]) {
                    used[graph::canonical_idx(e)] = 1;
                    dfs(g.edge(e).to);
                    trail.push_back(e);
                }
            }
        };
        dfs(v0);
        std::ranges::reverse(trail);
        return trail;
    }
    template<graph_type graph>
    std::optional<std::pair<node_index, std::vector<edge_index>>> euler_trail(graph const& g) {
        auto v0 = euler_start(g);
        if (!v0) {
            return std::nullopt;
        }
        auto result = try_euler_trail(g, *v0);
        if ((edge_index)result.size() != g.m()) {
            return std::nullopt;
        }
        return {{*v0, std::move(result)}};
    }
}
#endif // CP_ALGO_GRAPH_EULER_HPP
