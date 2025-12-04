#ifndef CP_ALGO_GRAPH_CYCLE_HPP
#define CP_ALGO_GRAPH_CYCLE_HPP
#include "dfs.hpp"
#include "base.hpp"
#include <deque>
namespace cp_algo::graph {
    template<graph_type graph>
    struct cycle_context: dfs_context<graph> {
        using base = dfs_context<graph>;
        using base::base;
        std::deque<edge_index> cycle;
        bool closed = false;
        int v0;
        
        void on_return_from_child(node_index v, edge_index e) {
            if (!empty(cycle) && !closed) {
                cycle.push_front(e);
                closed |= v == v0;
            }
        }
        
        void on_back_edge(node_index v, edge_index e) {
            if (empty(cycle)) {
                v0 = base::g->edge(e).to;
                base::done = true;
                closed = v == v0;
                cycle.push_front(e);
            }
        }
    };
    
    template<graph_type graph>
    std::deque<edge_index> find_cycle(graph const& g) {
        auto context = dfs<cycle_context>(g);
        return context.cycle;
    }
}
#endif // CP_ALGO_GRAPH_CYCLE_HPP
