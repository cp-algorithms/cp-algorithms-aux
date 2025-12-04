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
        
        auto const& adj = g.incidence_lists();
        
        struct frame {
            node_index v;
            edge_index ep;
            int sv; // edge index in stack_union
            enum { INIT, PROCESS_EDGES, HANDLE_CHILD } state;
        };
        
        std::stack<frame> dfs_stack;
        
        for (auto root: g.nodes()) {
            if (context.state[root] != unvisited) continue;
            
            dfs_stack.push({root, -1, 0, frame::INIT});
            
            while (!dfs_stack.empty()) {
                auto& f = dfs_stack.top();
                
                if (f.state == frame::INIT) {
                    context.state[f.v] = visiting;
                    context.tin[f.v] = context.low[f.v] = context.timer++;
                    context.stack.push(f.v);
                    f.sv = adj.head[f.v];
                    f.state = frame::PROCESS_EDGES;
                }
                
                if (f.state == frame::HANDLE_CHILD) {
                    auto e = adj.data[f.sv];
                    f.sv = adj.next[f.sv];
                    node_index u = g.edge(e).to;
                    context.low[f.v] = std::min(context.low[f.v], context.low[u]);
                    context.on_tree_edge(f.v, u);
                    f.state = frame::PROCESS_EDGES;
                }
                
                // PROCESS_EDGES
                bool found_child = false;
                while (f.sv != 0) {
                    auto e = adj.data[f.sv];
                    
                    if constexpr (undirected_graph_type<graph>) {
                        if (f.ep == graph::opposite_idx(e)) {
                            f.sv = adj.next[f.sv];
                            continue;
                        }
                    }
                    
                    node_index u = g.edge(e).to;
                    if (context.state[u] == unvisited) {
                        f.state = frame::HANDLE_CHILD;
                        dfs_stack.push({u, e, 0, frame::INIT});
                        found_child = true;
                        break;
                    } else if (context.state[u] != blocked) {
                        context.low[f.v] = std::min(context.low[f.v], context.tin[u]);
                    }
                    f.sv = adj.next[f.sv];
                }
                
                if (found_child) continue;
                
                // All edges processed
                context.state[f.v] = visited;
                context.on_exit(f.v);
                dfs_stack.pop();
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
