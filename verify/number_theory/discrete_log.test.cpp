// @brief Discrete Logarithm
#define PROBLEM "https://judge.yosupo.jp/problem/discrete_logarithm_mod"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/number_theory/discrete_log.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo;
using namespace math;
using base = dynamic_modint<>;

void solve() {
    int x, y, m;
    cin >> x >> y >> m;
    auto res = discrete_log(x, y, m);
    if(res) {
        cout << *res << "\n";
    } else {
        cout << -1 << "\n";
    }
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
