// @brief Cycle Detection (Directed)
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/cycle_detection
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("tune=native")
#include "cp-algo/graph/cycle.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, m;
    cin >> n >> m;
    graph<directed> g(n);
    g.read_edges(m);
    auto res = find_cycle(g);
    if(empty(res)) {
        cout << -1 << "\n";
    } else {
        cout << size(res) << "\n";
        for(auto it: res) {cout << it / 2 << '\n';}
    }
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    while(t--) {
        solve();
    }
}
