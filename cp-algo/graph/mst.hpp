#ifndef CP_ALGO_GRAPH_MST_HPP
#define CP_ALGO_GRAPH_MST_HPP
#include "base.hpp"
#include "../structures/dsu.hpp"
#include <algorithm>
namespace cp_algo::graph {
    template<weighted_edge_type edge_t>
    auto mst(graph<undirected, edge_t> const& g) {
        std::vector<std::pair<int64_t, edge_index>> edges;
        g.call_edges([&](edge_index e) {
            edges.emplace_back(g.edge(e).w, e);
        });
        std::ranges::sort(edges);
        structures::dsu me(g.n());
        int64_t total = 0;
        std::vector<edge_index> mst;
        for(auto [w, e]: edges) {
            if(me.uni(g.edge(e ^ 1).to, g.edge(e).to)) {
                total += w;
                mst.push_back(e);
            }
        }
        return std::pair{total, mst};
    }
}
#endif // CP_ALGO_GRAPH_MST_HPP
