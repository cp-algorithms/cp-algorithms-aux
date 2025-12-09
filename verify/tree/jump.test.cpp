// @brief Jump on Tree
#define PROBLEM "https://judge.yosupo.jp/problem/jump_on_tree"
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
    g.read_edges(n - 1);
    heavy_light hld(g);
    for(int i = 0; i < q; i++) {
        node_index u, v;
        int k;
        cin >> u >> v >> k;
        auto ores = hld.jump(u, v, k);
        if (!ores) {
            cout << -1 << '\n';
        } else {
            cout << *ores << '\n';
        }
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
