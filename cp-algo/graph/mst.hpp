#ifndef CP_ALGO_GRAPH_MST_HPP
#define CP_ALGO_GRAPH_MST_HPP
#include "base.hpp"
#include "../structures/dsu.hpp"
#include "../util/sort.hpp"
#include <algorithm>
namespace cp_algo::graph {
    template<weighted_undirected_graph_type graph>
    std::pair<int64_t, std::vector<edge_index>> mst(graph const& g) {
        struct edge {
            int64_t w;
            edge_index i;
        };
        auto edges = g.edge_indices() | std::ranges::to<std::vector>();
        radix_sort(edges, [&](auto e) {
            return g.edge(e).w;
        });
        structures::dsu me(g.n());
        int64_t total = 0;
        std::vector<edge_index> mst;
        for(auto e: edges) {
            if(me.uni(g.edge(e).to, g.edge(graph::opposite_idx(e)).to)) {
                total += g.edge(e).w;
                mst.push_back(graph::canonical_idx(e));
            }
        }
        return {total, mst};
    }
}
#endif // CP_ALGO_GRAPH_MST_HPP
