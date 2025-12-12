// @brief Minimum Spanning Tree
#define PROBLEM "https://judge.yosupo.jp/problem/minimum_spanning_tree"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/graph/mst.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, m;
    cin >> n >> m;
    weighted_graph g(n);
    g.read_edges(m);
    auto [X, E] = mst(g);
    cout << X << "\n";
    for(int e: E) {cout << e << " ";}
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
