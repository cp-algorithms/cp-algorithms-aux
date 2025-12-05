// @brief Eulerian Trail (Undirected)
#define PROBLEM "https://judge.yosupo.jp/problem/eulerian_trail_undirected"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/graph/euler.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, m;
    cin >> n >> m;
    graph g(n);
    g.read_edges(m);
    auto trail = euler_trail(g);
    if(!trail) {
        cout << "No" << "\n";
    } else {
        auto [v, es] = *trail;
        cout << "Yes" << "\n";
        cout << v;
        for(auto e: es) {
            v = g.edge(e).traverse(v);
            cout << ' ' << v;
        }
        cout << "\n";
        for(auto e: es) {cout << e << ' ';}
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
