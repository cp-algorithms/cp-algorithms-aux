#ifndef CP_ALGO_GRAPH_XOR_DFS_HPP
#define CP_ALGO_GRAPH_XOR_DFS_HPP
#include "base.hpp"
#include <vector>

namespace cp_algo::graph {

/**
 * @brief XOR DFS on trees - processes leaves to root using XOR trick
 * 
 * Uses XOR to track remaining edges instead of storing parent pointers.
 * Calls callback for each leaf being processed, passing the leaf node and its parent edge.
 * 
 * @param tree The input tree (undirected graph)
 * @param callback Function called for each leaf with (node, edge) parameters
 * @param root Root node (default 0)
 * @return Vector where result[v] = XOR of all incident edge indices at node v
 * 
 * Time complexity: O(n)
 * Space complexity: O(n)
 */
template<undirected_graph_type graph>
auto xor_dfs(graph const& tree, auto &&callback, node_index root = 0) {
    std::vector<int> degree(tree.n());
    std::vector<edge_index> neig_xor(tree.n());
    for(auto v: tree.nodes()) {
        for(auto e: tree.outgoing(v)) {
            degree[v]++;
            neig_xor[v] ^= e;
        }
    }
    neig_xor[root] ^= edge_index(-1);
    degree[root]++;
    for(auto v: tree.nodes()) {
        while(degree[v] == 1) {
            edge_index ep = neig_xor[v];
            callback(v, ep);
            degree[v]--;
            if (v == root) break;
            v = tree.edge(ep).traverse(v);
            degree[v]--;
            neig_xor[v] ^= ep;
        }
    }
    return neig_xor;
}

} // namespace cp_algo::graph
#endif // CP_ALGO_GRAPH_XOR_DFS_HPP
