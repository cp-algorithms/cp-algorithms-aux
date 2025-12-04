// @brief Cycle Detection (Undirected)
#define PROBLEM "https://judge.yosupo.jp/problem/cycle_detection_undirected"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/graph/cycle.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::graph;

void solve() {
    int n, m;
    cin >> n >> m;
    graph g(n);
    g.read_edges(m);
    auto res = find_cycle(g);
    if(empty(res)) {
        cout << -1 << "\n";
    } else {
        cout << size(res) << "\n";
        ranges::rotate(res, prev(end(res)));
        for(auto it: res) {cout << g.edge(it).to << ' ';}
        cout << "\n";
        ranges::rotate(res, next(begin(res)));
        for(auto it: res) {cout << graph<>::canonical_idx(it) << ' ';}
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
