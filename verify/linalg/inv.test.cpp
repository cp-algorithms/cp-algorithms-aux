// @brief Inverse Matrix
#define PROBLEM "https://judge.yosupo.jp/problem/inverse_matrix"
#pragma GCC optimize("Ofast,unroll-loops")
#include <bits/stdc++.h>
#include "cp-algo/linalg/matrix.hpp"

using namespace std;
using namespace cp_algo::linalg;
using namespace cp_algo::math;

const int64_t mod = 998244353;

void solve() {
    int n;
    cin >> n;
    matrix<modint<mod>> a(n, n);
    a.read();
    auto [d, ai] = a.inv();
    if(d == 0) {
        cout << -1 << "\n";
    } else {
        ai.print();
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
