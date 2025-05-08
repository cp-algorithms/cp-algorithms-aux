// @brief Minimum Spanning Tree
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/minimum_spanning_tree
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/graph/mst.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, m;
    cin >> n >> m;
    graph<undirected, weighted_edge> g(n);
    g.read_edges(m);
    auto [X, E] = mst(g);
    cout << X << "\n";
    for(int e: E) {cout << e / 2 << " ";}
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
