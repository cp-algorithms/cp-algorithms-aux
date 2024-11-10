// @brief Counting Spanning Trees (Directed)
#define PROBLEM "https://judge.yosupo.jp/problem/counting_spanning_tree_directed"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("tune=native")
#include "cp-algo/linalg/matrix.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;
using namespace cp_algo::linalg;

const int64_t mod = 998244353;
using base = modint<mod>;

void solve() {
    int n, m, r;
    cin >> n >> m >> r;
    matrix<base> a(n);
    for(int i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        a[u][v] -= 1;
        a[v][v] += 1;
    }
    for(int i = 0; i < n; i++) {
        a[r][i] = a[i][r] = 0;
    }
    a[r][r] = 1;
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
