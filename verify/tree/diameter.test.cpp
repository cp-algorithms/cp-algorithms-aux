// @brief Tree Diameter
#define PROBLEM "https://judge.yosupo.jp/problem/tree_diameter"
#pragma GCC optimize("Ofast,unroll-loops")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/tree/diameter.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n;
    cin >> n;
    weighted_graph g(n);
    g.read_edges(n - 1);
    auto [d, v, path] = tree_diameter<diameter_mode::recover_path>(g);
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
