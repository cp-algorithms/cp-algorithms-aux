// @brief Lowest Common Ancestor
#define PROBLEM "https://judge.yosupo.jp/problem/lca"
#pragma GCC optimize("O3,unroll-loops")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/tree/hld.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, q;
    cin >> n >> q;
    graph g(n);
    vector<edge_index> p(n);
    p[0] = -1;
    for(int i = 1; i < n; i++) {
        node_index v;
        cin >> v;
        p[i] = g.add_edge(v, i);
    }
    heavy_light hld(g, 0, &p);
    for(int i = 0; i < q; i++) {
        node_index u, v;
        cin >> u >> v;
        cout << hld.lca(u, v) << '\n';
    }
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
