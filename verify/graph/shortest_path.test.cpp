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
    weighted_digraph g(n);
    g.read_edges(m);
    auto opath = shortest_path(g, s, t);
    if(!opath) {
        cout << -1 << "\n";
    } else {
        auto [d, path] = *opath;
        cout << d << ' ' << size(path) << "\n";
        int v = s;
        for(auto e: path) {
            int w = g.edge(e).traverse(v);
            cout << exchange(v, w) << ' ' << w << "\n";
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
