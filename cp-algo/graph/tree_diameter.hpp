#ifndef CP_ALGO_GRAPH_TREE_DIAMETER_HPP
#define CP_ALGO_GRAPH_TREE_DIAMETER_HPP
#include "xor_dfs.hpp"
#include "base.hpp"
#include <tuple>
#include <string>
#include <algorithm>

namespace cp_algo::graph {

/**
 * @brief Finds the diameter of a weighted tree
 * 
 * Returns the maximum distance between any two nodes in the tree.
 * Uses XOR DFS to traverse the tree while tracking path distances.
 * 
 * @param g The input weighted tree
 * @return Tuple of (diameter_length, starting_node, path_edges)
 *         The path is represented as a sequence of edge indices from start to end
 * 
 * Time complexity: O(n)
 * Space complexity: O(n)
 */
template<weighted_undirected_graph_type graph>
std::tuple<int64_t, node_index, std::basic_string<edge_index>> tree_diameter(graph const& g) {
    struct up_path {
        int64_t length = 0;
        node_index start;
    };
    std::vector<up_path> up(g.n());
    for(auto v: g.nodes()) {
        up[v].start = v;
    }
    up_path s, t;
    auto parents = xor_dfs(g, [&](node_index v, edge_index ep) {
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
    });
    auto collect = [&](up_path v) {
        std::basic_string<edge_index> path;
        while(v.length) {
            edge_index ep = parents[v.start];
            path.push_back(ep);
            v.length -= g.edge(ep).w;
            v.start = g.edge(ep).traverse(v.start);
        }
        return path;
    };
    auto paths = collect(s);
    auto patht = collect(t);
    std::ranges::reverse(patht);
    return {s.length + t.length, s.start, paths += patht};
}

} // namespace cp_algo::graph
#endif // CP_ALGO_GRAPH_TREE_DIAMETER_HPP
