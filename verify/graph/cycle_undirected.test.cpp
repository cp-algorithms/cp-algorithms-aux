// @brief Cycle Detection (Undirected)
#define PROBLEM "https://judge.yosupo.jp/problem/cycle_detection_undirected"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/graph/cycle.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, m;
    cin >> n >> m;
    graph g(n);
    g.read_edges(m);
    auto [v, cycle] = find_cycle(g);
    if(empty(cycle)) {
        cout << -1 << "\n";
    } else {
        cout << size(cycle) << "\n";
        cout << v;
        for(auto it: cycle | views::take(size(cycle) - 1)) {
            v = g.edge(it).traverse(v);
            cout << ' ' << v;
        }
        cout << "\n";
        for(auto it: cycle) {cout << it << ' ';}
        cout << "\n";
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
