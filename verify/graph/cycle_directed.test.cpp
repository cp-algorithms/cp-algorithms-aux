// @brief Cycle Detection (Directed)
#define PROBLEM "https://judge.yosupo.jp/problem/cycle_detection"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/graph/cycle.hpp"

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, m;
    cin >> n >> m;
    digraph g(n);
    g.read_edges(m);
    auto [v, cycle] = find_cycle(g);
    if(empty(cycle)) {
        cout << -1 << "\n";
    } else {
        cout << size(cycle) << "\n";
        for(auto it: cycle) {cout << it << '\n';}
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
