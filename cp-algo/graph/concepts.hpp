#ifndef CP_ALGO_GRAPH_CONCEPTS_HPP
#define CP_ALGO_GRAPH_CONCEPTS_HPP
#include "edge_types.hpp"
#include <type_traits>

namespace cp_algo::graph {
    // Shared graph mode enum for all graph headers
    enum graph_mode { directed, undirected };
    // Traits: true for types that expose `edge_t` and static `mode`
    template<typename T, typename = void>
    struct graph_traits : std::false_type {};

    template<typename T>
    struct graph_traits<T, std::void_t<typename T::edge_t, decltype(T::mode)>> : std::true_type {
        using edge_t = typename T::edge_t;
        static constexpr auto mode = T::mode;
        static constexpr bool is_directed = mode == directed;
        static constexpr bool is_undirected = mode == undirected;
        static constexpr bool is_weighted = weighted_edge_type<edge_t>;
    };

    // Concepts
    template<typename G>
    concept graph_type = graph_traits<G>::value;

    template<typename G>
    concept digraph_type = graph_type<G> && graph_traits<G>::is_directed;

    template<typename G>
    concept undirected_graph_type = graph_type<G> && graph_traits<G>::is_undirected;

    template<typename G>
    concept weighted_graph_type = graph_type<G> && graph_traits<G>::is_weighted;

    template<typename G>
    concept weighted_digraph_type = digraph_type<G> && graph_traits<G>::is_weighted;

    template<typename G>
    concept weighted_undirected_graph_type = undirected_graph_type<G> && graph_traits<G>::is_weighted;
}
#endif // CP_ALGO_GRAPH_CONCEPTS_HPP
