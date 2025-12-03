// @brief Biconnected Components
#define PROBLEM "https://judge.yosupo.jp/problem/biconnected_components"
#pragma GCC optimize("Ofast,unroll-loops")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/graph/bcc.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, m;
    cin >> n >> m;
    graph g(n);
    g.read_edges(m);
    auto comps = biconnected_components(g);
    cout << size(comps) << '\n';
    for(auto const& comp: comps.rows()) {
        cout << size(comp) << ' ';
        for(auto v: comp) {
            cout << v << ' ';
        }
        cout << '\n';
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
