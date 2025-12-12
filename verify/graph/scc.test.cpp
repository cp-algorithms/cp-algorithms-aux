// @brief Strongly Connected Components
#define PROBLEM "https://judge.yosupo.jp/problem/scc"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <iostream>
//#include "blazingio/blazingio.min.hpp"
#include "cp-algo/graph/tarjan.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, m;
    cin >> n >> m;
    digraph g(n);
    g.read_edges(m);
    auto comps = strongly_connected_components(g);
    cout << size(comps) << '\n';
    for(auto const& comp: comps.rows() | views::reverse) {
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
