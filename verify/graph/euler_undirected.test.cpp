// @brief Eulerian Trail (Undirected)
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/eulerian_trail_undirected
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("tune=native")
#include "cp-algo/graph/euler.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, m;
    cin >> n >> m;
    graph<undirected> g(n);
    g.read_edges(m);
    auto [v0, es] = euler_trail(g);
    if(ssize(es) != m) {
        cout << "No" << "\n";
    } else {
        cout << "Yes" << "\n";
        cout << v0 << ' ';
        for(auto e: es) {cout << g.edge(e).to << ' ';}
        cout << "\n";
        for(auto e: es) {cout << e / 2 << ' ';}
        cout << "\n";
    }
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    cin >> t;
    while(t--) {
        solve();
    }
}
