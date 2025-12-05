// @brief Minimum Diameter Spanning Tree
#define PROBLEM "https://judge.yosupo.jp/problem/minimum_diameter_spanning_tree"
#pragma GCC optimize("Ofast,unroll-loops")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/graph/shortest_path.hpp"
#include "cp-algo/tree/diameter.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

template<weighted_undirected_graph_type graph>
auto minimum_diameter_spanning_tree(graph const& g) {
    int64_t min_diameter = shortest_path_context::inf;
    std::vector<edge_index> best_edges;
    for(auto v: g.nodes()) {
        auto [_, p] = dijkstra(g, v);
        auto d = tree_diameter(g, &p);
        if (d < min_diameter) {
            min_diameter = d;
            best_edges = p;
        }
    }
    return std::pair{min_diameter, best_edges};
}

void solve() {
    int n, m;
    cin >> n >> m;
    weighted_graph g(n);
    g.read_edges(m);
    auto [X, E] = minimum_diameter_spanning_tree(g);
    cout << X << '\n';
    for(auto e: E) {
        if (e != -1) {
            cout << e << " ";
        }
    }
    cout << "\n";
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    //cin >> t;
    while(t--) {
        solve();
    }
}
