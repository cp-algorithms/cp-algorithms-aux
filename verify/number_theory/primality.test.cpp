// @brief Primality Test
#define PROBLEM "https://judge.yosupo.jp/problem/primality_test"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/math/primality.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

void solve() {
    int64_t m;
    cin >> m;
    cout << (is_prime(m) ? "Yes" : "No") << "\n";
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
