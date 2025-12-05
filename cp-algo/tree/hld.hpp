#ifndef CP_ALGO_TREE_HLD_HPP
#define CP_ALGO_TREE_HLD_HPP
#include "ascending_dfs.hpp"
#include "../graph/base.hpp"
#include "../util/big_alloc.hpp"
#include <vector>
#include <ranges>

namespace cp_algo::graph {
    struct heavy_light {
        big_vector<node_index> size, in, up, par;
        
        template<undirected_graph_type graph>
        heavy_light(graph const& g, node_index root = 0, std::vector<edge_index> const* parents_ptr = nullptr):
            size(g.n(), 1), in(g.n()), up(g.n()), par(g.n()) {
            big_vector<node_index> topsort;
            topsort.reserve(g.n());
            auto push_size = [&](node_index v, edge_index e) {
                topsort.push_back(v);
                if (v != root) {
                    auto p = g.edge(e).traverse(v);
                    size[p] += size[v];
                    par[v] = p;
                }
            };
            if (parents_ptr) {
                parent_dfs(g, *parents_ptr, push_size);
            } else {
                xor_dfs(g, push_size, root);
            }
            par[root] = up[root] = root;
            for(auto v: topsort | std::views::reverse) {
                if (size[v] == 1) continue;
                node_index big = -1;
                for(auto e: g.outgoing(v)) {
                    auto u = g.edge(e).traverse(v);
                    if (size[u] > size[v]) continue;
                    if (big == -1 || size[u] > size[big]) {
                        big = u;
                    }
                }
                int t = in[v] + size[big];
                for(auto e: g.outgoing(v)) {
                    auto u = g.edge(e).traverse(v);
                    if (size[u] > size[v]) continue;
                    if (u == big) {
                        in[u] = in[v] + 1;
                        up[u] = up[v];
                    } else {
                        in[u] = t + 1;
                        t += size[u];
                        up[u] = u;
                    }
                }
            }
        }
        node_index lca(node_index a, node_index b) {
            while (up[a] != up[b]) {
                if (in[up[a]] < in[up[b]]) {
                    b = par[up[b]];
                } else {
                    a = par[up[a]];
                }
            }
            return in[a] < in[b] ? a : b;
        }
    };
}
#endif // CP_ALGO_TREE_HLD_HPP
