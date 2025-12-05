#ifndef CP_ALGO_TREE_HLD_HPP
#define CP_ALGO_TREE_HLD_HPP
#include "ascending_dfs.hpp"
#include "../graph/base.hpp"
#include "../util/big_alloc.hpp"
namespace cp_algo::graph {
    struct heavy_light {
        big_vector<node_index> size, in, up, par;
        
        template<undirected_graph_type graph>
        heavy_light(graph const& g, node_index root = 0, std::vector<edge_index> const* parents_ptr = nullptr):
            size(g.n(), 1), in(g.n()), up(g.n()), par(g.n()) {
            big_vector<node_index> topsort;
            topsort.reserve(g.n());
            auto push_size = [&](node_index v, edge_index e) {
                if (size[v] > 1) {
                    topsort.push_back(v);
                }
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
        enum lca_mode { without_distances, with_distances };
        template<lca_mode mode = without_distances>
        auto lca(node_index a, node_index b) {
            int dista = 0, distb = 0;
            while (up[a] != up[b]) {
                if (in[up[a]] < in[up[b]]) {
                    if constexpr (mode == with_distances) distb += in[b] - in[up[b]] + 1;
                    b = par[up[b]];
                } else {
                    if constexpr (mode == with_distances) dista += in[a] - in[up[a]] + 1;
                    a = par[up[a]];
                }
            }
            node_index c = in[a] < in[b] ? a : b;
            if constexpr (mode == with_distances) {
                return std::tuple{c, dista + in[a] - in[c], distb + in[b] - in[c]};
            } else {
                return c;
            }
        }
        big_vector<node_index> rin;
        void compute_rin() {
            if (empty(rin)) {
                rin.resize(std::size(in));
                for (auto [v, inv]: in | std::views::enumerate) {
                    rin[inv] = node_index(v);
                }
            }
        }
        node_index jump_up(node_index v, int steps) {
            compute_rin();
            while (steps > 0) {
                int path_dist = in[v] - in[up[v]];
                if (steps <= path_dist) {
                    return rin[in[v] - steps];
                }
                steps -= path_dist + 1;
                v = par[up[v]];
            }
            return v;
        }
        std::optional<node_index> jump(node_index from, node_index to, int steps) {
            compute_rin();
            auto [l, dist_from, dist_to] = lca<with_distances>(from, to);
            auto dist = dist_from + dist_to;
            if (steps > dist) return std::nullopt;
            if (steps <= dist_from) {
                return jump_up(from, steps);
            } else {
                return jump_up(to, dist - steps);
            }
        }
    };
}
#endif // CP_ALGO_TREE_HLD_HPP
