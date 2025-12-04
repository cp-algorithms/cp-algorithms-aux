#ifndef CP_ALGO_GRAPH_EULER_HPP
#define CP_ALGO_GRAPH_EULER_HPP
#include "base.hpp"
#include <algorithm>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>
#include <stack>
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
    std::deque<edge_index> try_euler_trail(graph const& g, node_index v0) {
        std::deque<edge_index> trail;
        enum state { unvisited, visited };
        std::vector<state> state(g.m());
        auto const& adj = g.incidence_lists();
        auto head = adj.head;
        
        std::stack<edge_index> stack;
        stack.push(-1);
        
        while (!empty(stack)) {
            auto ep = stack.top();
            node_index w = ~ep ? g.edge(ep).to : v0;
            bool found_edge = false;
            
            while (head[w] != 0) {
                auto e = adj.data[std::exchange(head[w], adj.next[head[w]])];
                if(state[graph::canonical_idx(e)] == unvisited) {
                    state[graph::canonical_idx(e)] = visited;
                    stack.push(e);
                    found_edge = true;
                    break;
                }
            }
            
            if (!found_edge) {
                stack.pop();
                if (~ep) {
                    trail.push_front(ep);
                }
            }
        }
        return trail;
    }
    template<graph_type graph>
    std::optional<std::pair<node_index, std::deque<edge_index>>> euler_trail(graph const& g) {
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
