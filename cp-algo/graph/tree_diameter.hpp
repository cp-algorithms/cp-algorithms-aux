#ifndef CP_ALGO_GRAPH_TREE_DIAMETER_HPP
#define CP_ALGO_GRAPH_TREE_DIAMETER_HPP
#include "ascending_dfs.hpp"
#include "base.hpp"
#include <tuple>
#include <string>
#include <algorithm>
#include <memory>

namespace cp_algo::graph {
    enum diameter_mode { recover_path, no_recover };
    template<diameter_mode mode = no_recover, weighted_undirected_graph_type graph>
    auto tree_diameter(graph const& g, std::shared_ptr<const std::vector<edge_index>> parents = nullptr) {
        struct up_path {
            int64_t length = 0;
            node_index start = 0;
        };
        std::vector<up_path> up(g.n());
        for(auto v: g.nodes()) {
            up[v].start = v;
        }
        up_path s, t;
        auto callback = [&](node_index v, edge_index ep) {
            if (ep == edge_index(-1)) return;
            node_index u = g.edge(ep).traverse(v);
            up[v].length += g.edge(ep).w;
            if (up[v].length + up[u].length > s.length + t.length) {
                s = up[v];
                t = up[u];
            }
            if (up[v].length > up[u].length) {
                up[u] = up[v];
            }
        };
        if (parents) {
            parent_dfs(g, *parents, callback);
        } else {
            parents = std::make_shared<std::vector<edge_index>>(xor_dfs(g, callback));
        }
        if constexpr (mode == no_recover) {
            return s.length + t.length;
        }
        auto collect = [&](up_path v) {
            std::basic_string<edge_index> path;
            while(v.length) {
                edge_index ep = (*parents)[v.start];
                path.push_back(ep);
                v.length -= g.edge(ep).w;
                v.start = g.edge(ep).traverse(v.start);
            }
            return path;
        };
        auto paths = collect(s);
        auto patht = collect(t);
        std::ranges::reverse(patht);
        return std::tuple{s.length + t.length, s.start, paths += patht};
    }
}
#endif // CP_ALGO_GRAPH_TREE_DIAMETER_HPP
