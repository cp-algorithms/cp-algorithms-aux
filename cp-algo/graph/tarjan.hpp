#ifndef CP_ALGO_GRAPH_TARJAN_HPP
#define CP_ALGO_GRAPH_TARJAN_HPP
#include "base.hpp"
#include "../structures/csr.hpp"
#include <algorithm>
#include <cassert>
#include <stack>
namespace cp_algo::graph {
    enum node_state { unvisited, visiting, visited, blocked };
    template<graph_type graph>
    struct tarjan_context {
        big_vector<int> tin, low;
        big_vector<node_state> state;
        std::stack<int> stack;
        graph const* g;
        int timer;
        structures::csr<node_index> components;
        tarjan_context(graph const& g):
            tin(g.n()), low(g.n()), state(g.n()), g(&g), timer(0) {
            components.reserve_data(g.n());
        }

        void on_tree_edge(node_index, edge_index) {}
        void on_exit(node_index) {}

        void collect(node_index v) {
            components.new_row();
            node_index u;
            do {
                u = stack.top();
                stack.pop();
                state[u] = blocked;
                components.push(u);
            } while(u != v);
        }
    };
    template<template<typename> class Context, graph_type graph>
    auto tarjan(graph const& g) {
        Context<graph> context(g);
        auto dfs = [&](this auto &&dfs, node_index v, edge_index ep = -1) -> void {
            context.state[v] = visiting;
            context.tin[v] = context.low[v] = context.timer++;
            context.stack.push(v);
            for(auto e: g.outgoing(v)) {
                if constexpr (undirected_graph_type<graph>) {
                    if (ep == graph::opposite_idx(e)) {
                        continue;
                    }
                }
                node_index u = g.edge(e).to;
                if (context.state[u] == unvisited) {
                    dfs(u, e);
                    context.low[v] = std::min(context.low[v], context.low[u]);
                    context.on_tree_edge(v, u);
                } else if (context.state[u] != blocked) {
                    context.low[v] = std::min(context.low[v], context.tin[u]);
                }
            }
            context.state[v] = visited;
            context.on_exit(v);
        };
        for (auto v: g.nodes()) {
            if (context.state[v] == unvisited) {
                dfs(v);
            }
        }
        return context.components;
    }
    template<graph_type graph>
    struct exit_context: tarjan_context<graph> {
        using base = tarjan_context<graph>;
        using base::base;
        
        void on_exit(node_index v) {
            if (base::low[v] == base::tin[v]) {
                base::collect(v);
            }
        }
    };
    // Tarjan's algorithm for Strongly Connected Components
    // returns components in reverse topological order
    template<digraph_type graph>
    auto strongly_connected_components(graph const& g) {
        return tarjan<exit_context>(g);
    }
    
    // Tarjan's algorithm for Two-Edge-Connected Components
    template<undirected_graph_type graph>
    auto two_edge_connected_components(graph const& g) {
        return tarjan<exit_context>(g);
    }

    template<undirected_graph_type graph>
    struct bcc_context: tarjan_context<graph> {
        using base = tarjan_context<graph>;
        using base::base;
        void on_tree_edge(node_index v, node_index u) {
            if (base::low[u] >= base::tin[v]) {
                base::collect(u);
                base::components.push(v);
            }
        }
        void on_exit(node_index v) {
            if (std::empty(base::g->outgoing(v))) {
                base::collect(v);
            }
        }
    };
    // Tarjan's algorithm for Biconnected Components (vertex-biconnected)
    template<undirected_graph_type graph>
    auto biconnected_components(graph const& g) {
        return tarjan<bcc_context>(g);
    }
}
#endif // CP_ALGO_GRAPH_TARJAN_HPP
