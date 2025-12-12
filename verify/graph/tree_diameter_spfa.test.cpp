// @brief Tree Diameter (SPFA)
#define PROBLEM "https://judge.yosupo.jp/problem/tree_diameter"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/graph/shortest_path.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

template<weighted_undirected_graph_type graph>
std::tuple<int64_t, node_index, cp_algo::big_vector<edge_index>> tree_diameter(graph const& g) {
    auto [d1, _] = spfa(g, 0);
    node_index s = 0;
    for(auto v: g.nodes()) {
        if (d1[v] > d1[s]) {
            s = v;
        }
    }
    auto [d2, pre] = spfa(g, s);
    node_index t = 0;
    for(auto v: g.nodes()) {
        if (d2[v] > d2[t]) {
            t = v;
        }
    }
    return {d2[t], s, recover_path(g, pre, s, t)};
}

void solve() {
    int n;
    cin >> n;
    weighted_graph g(n);
    g.read_edges(n - 1);
    auto [d, v, path] = tree_diameter(g);
    cout << d << ' ' << size(path) + 1 << '\n';
    cout << v;
    for(auto e: path) {
        v = g.edge(e).traverse(v);
        cout << ' ' << v;
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
