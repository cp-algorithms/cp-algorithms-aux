#ifndef CP_ALGO_GRAPH_BASE_HPP
#define CP_ALGO_GRAPH_BASE_HPP
#include "edge_types.hpp"
#include "concepts.hpp"
#include "../structures/stack_union.hpp"
#include <ranges>
namespace cp_algo::graph {
    using edge_index = int;
    template<edge_type _edge_t = edge_base, graph_mode _mode = undirected>
    struct graph {
        using edge_t = _edge_t;
        static constexpr auto mode = _mode;
        using incidence_list = structures::stack_union<edge_index>;
        graph(int n, int v0 = 0): v0(v0), adj(n) {}

        void add_edge(node_index u, edge_t e) {
            adj.push(u, (edge_index)size(E));
            E.push_back(e);
            if constexpr (mode == undirected) {
                adj.push(e.to, (edge_index)size(E));
                E.push_back(e.backedge(u));
            }
        }
        void add_edge(node_index u, auto... Args) {
            add_edge(u, edge_t(Args...));
        }
        void read_edges(node_index m) {
            adj.reserve(mode == undirected ? 2 * m : m);
            for(edge_index i = 0; i < m; i++) {
                auto [u, e] = edge_t::read(v0);
                add_edge(u, e);
            }
        }
        auto outgoing(node_index v) const {
            return adj[v];
        }
        auto nodes() const {
            return std::views::iota(node_index(0), n());
        }

        auto edges() const {
            if constexpr (mode == undirected) {
                return E | std::views::stride(2);
            } else {
                return E | std::views::all;
            }
        }
        auto edge_indices() const {
            auto indices = std::views::iota(edge_index(0), edge_index(size(E)));
            if constexpr (mode == undirected) {
                return indices | std::views::stride(2);
            } else {
                return indices;
            }
        }
        incidence_list const& incidence_lists() const {return adj;}
        edge_t const& edge(edge_index e) const {return E[e];}
        node_index n() const {return (node_index)adj.size();}
        edge_index m() const {
            return (edge_index)size(edges());
        }
        static edge_index canonical_idx(edge_index e) {
            if constexpr (mode == undirected) {
                return e / 2;
            } else {
                return e;
            }
        }
        static edge_index opposite_idx(edge_index e) {
            static_assert(mode == undirected, "opposite_index is only defined for undirected graphs");
            return e ^ 1;
        }
    private:
        node_index v0;
        big_vector<edge_t> E;
        incidence_list adj;
    };
    // aliases for most standard cases
    template<edge_type edge_t = edge_base>
    using digraph = graph<edge_t, directed>;
    template<weighted_edge_type edge_t = weighted_edge, graph_mode mode = undirected>
    using weighted_graph = graph<edge_t, mode>;
    template<weighted_edge_type edge_t = weighted_edge>
    using weighted_digraph = digraph<edge_t>;
}
#endif // CP_ALGO_GRAPH_BASE_HPP
