#ifndef CP_ALGO_GRAPH_TARJAN_HPP
#define CP_ALGO_GRAPH_TARJAN_HPP
#include "dfs.hpp"
#include "base.hpp"
#include "../structures/csr.hpp"
#include <algorithm>
#include <cassert>
#include <stack>
namespace cp_algo::graph {
    template<graph_type graph>
    struct tarjan_context: dfs_context<graph> {
        using base = dfs_context<graph>;
        big_vector<int> tin, low;
        std::stack<int> stack;
        int timer;
        structures::csr<node_index> components;
        
        tarjan_context(graph const& g): base(g),
            tin(g.n()), low(g.n()), timer(0) {
            components.reserve_data(g.n());
        }

        void on_enter(node_index v) {
            tin[v] = low[v] = timer++;
            stack.push(v);
        }

        void on_return_from_child(node_index v, edge_index e) {
            node_index u = base::g->edge(e).traverse(v);
            low[v] = std::min(low[v], low[u]);
        }

        void on_back_edge(node_index v, edge_index e) {
            node_index u = base::g->edge(e).traverse(v);
            low[v] = std::min(low[v], tin[u]);
        }

        void on_forward_cross_edge(node_index v, edge_index e) {
            node_index u = base::g->edge(e).traverse(v);
            low[v] = std::min(low[v], tin[u]);
        }

        void on_tree_edge_processed(node_index, node_index, edge_index) {}

        void collect(node_index v) {
            components.new_row();
            node_index u;
            do {
                u = stack.top();
                stack.pop();
                base::state[u] = blocked;
                components.push(u);
            } while(u != v);
        }
    };
    template<graph_type graph>
    struct exit_context: tarjan_context<graph> {
        using tarjan_context<graph>::tarjan_context;
        
        void on_exit(node_index v) {
            if (this->low[v] == this->tin[v]) {
                this->collect(v);
            }
        }
    };
    // Tarjan's algorithm for Strongly Connected Components
    // returns components in reverse topological order
    template<digraph_type graph>
    auto strongly_connected_components(graph const& g) {
        return dfs<exit_context>(g).components;
    }
    
    // Tarjan's algorithm for Two-Edge-Connected Components
    template<undirected_graph_type graph>
    auto two_edge_connected_components(graph const& g) {
        return dfs<exit_context>(g).components;
    }

    template<undirected_graph_type graph>
    struct bcc_context: tarjan_context<graph> {
        using base = tarjan_context<graph>;
        using base::base;
        
        void on_return_from_child(node_index v, edge_index e) {
            base::on_return_from_child(v, e);
            node_index u = base::g->edge(e).traverse(v);
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
        return dfs<bcc_context>(g).components;
    }
}
#endif // CP_ALGO_GRAPH_TARJAN_HPP
