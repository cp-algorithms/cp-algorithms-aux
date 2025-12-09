// @brief Sqrt Mod
#define PROBLEM "https://judge.yosupo.jp/problem/sqrt_mod"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("tune=native")
#include "cp-algo/number_theory/discrete_sqrt.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;
using base = dynamic_modint<>;

void solve() {
    int y, p;
    cin >> y >> p;
    base::switch_mod(p);
    auto res = sqrt(base(y));
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
