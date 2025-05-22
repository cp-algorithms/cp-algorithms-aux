// @brief Shortest Path
#define PROBLEM "https://judge.yosupo.jp/problem/shortest_path"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/graph/shortest_path.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, m, s, t;
    cin >> n >> m >> s >> t;
    graph<directed, weighted_edge> g(n);
    g.read_edges(m);
    auto [X, Y] = shortest_path(g, s, t);
    if(X == -1) {
        cout << X << "\n";
    } else {
        cout << X << ' ' << size(Y) << "\n";
        for(auto [u, e]: Y) {
            cout << u << ' ' << g.edge(e).to << "\n";
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
