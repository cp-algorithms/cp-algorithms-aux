// @brief Counting Spanning Trees (Undirected)
#define PROBLEM "https://judge.yosupo.jp/problem/counting_spanning_tree_undirected"
#pragma GCC optimize("Ofast,unroll-loops")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/combinatorics.hpp"
#include "cp-algo/linalg/matrix.hpp"

using namespace std;
using namespace cp_algo::math;
using namespace cp_algo::linalg;

const int64_t mod = 998244353;
using base = modint<mod>;

void solve() {
    int n, m;
    cin >> n >> m;
    matrix<base> a(n);
    for(int i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        a[u][v] -= 1;
        a[v][u] -= 1;
        a[v][v] += 1;
        a[u][u] += 1;
    }
    for(int i = 0; i < n; i++) {
        a[0][i] = a[i][0] = 0;
    }
    a[0][0] = 1;
    cout << a.det() << "\n";
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    // cin >> t;
    while(t--) {
        solve();
    }
}
