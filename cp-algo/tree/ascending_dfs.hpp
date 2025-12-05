#ifndef CP_ALGO_TREE_ASCENDING_DFS_HPP
#define CP_ALGO_TREE_ASCENDING_DFS_HPP

#include "../graph/base.hpp"
#include <cassert>
#include <vector>
#include <ranges>

namespace cp_algo::graph {
    // Generic ascending DFS that repeatedly peels leaves (degree 1).
    // `next(v)` must return the unique incident edge of `v` when degree[v] == 1.
    template<undirected_graph_type graph>
    void ascending_dfs(graph const& tree, auto &degree, auto &&next, auto &&callback, node_index root) {
        for (auto v : tree.nodes()) {
            while (degree[v] == 1) {
                edge_index ep = next(v);
                callback(v, ep);
                degree[v]--;
                v = tree.edge(ep).traverse(v);
                degree[v]--;
            }
        }
        callback(root, -1);
    }
    // XOR-based DFS that tracks the remaining incident edge via XOR.
    template<undirected_graph_type graph>
    auto xor_dfs(graph const& tree, auto &&callback, node_index root = 0) {
        std::vector<edge_index> neig_xor(tree.n());
        std::vector<int> degree(tree.n());
        for (auto v : tree.nodes()) {
            degree[v] = (int)std::ranges::distance(tree.outgoing(v));
        }
        degree[root] = 0;
        for (auto v : tree.nodes()) {
            for (auto e : tree.outgoing(v)) {
                neig_xor[v] ^= e;
            }
        }
        neig_xor[root] ^= edge_index(-1);
        ascending_dfs(tree, degree, [&](auto v) {
            edge_index ep = neig_xor[v];
            neig_xor[tree.edge(ep).traverse(v)] ^= ep;
            return ep;
        }, callback, root);

        return neig_xor; // parent edge for each node (root holds -1)
    }
    // DFS that uses a precomputed parent-edge array.
    template<undirected_graph_type graph>
    void parent_dfs(graph const& tree, std::vector<edge_index> const& parent, auto &&callback) {
        std::vector<int> degree(tree.n());
        node_index root = -1;
        for (auto [v, e] : parent | std::views::enumerate) {
            if (e != -1) {
                degree[v]++;
                degree[tree.edge(e).traverse(node_index(v))]++;
            } else {
                root = node_index(v);
            }
        }
        assert(root != -1);
        degree[root] = 0;
        ascending_dfs(tree, degree, [&](auto v) {
            return parent[v];
        }, callback, root);
    }
}
#endif // CP_ALGO_TREE_ASCENDING_DFS_HPP
