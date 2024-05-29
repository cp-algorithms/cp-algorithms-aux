// @brief Primitive Root
#define PROBLEM "https://judge.yosupo.jp/problem/primitive_root"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("tune=native")
#include "cp-algo/math/number_theory.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo;
using namespace math;
using base = dynamic_modint;

void solve() {
    int64_t p;
    cin >> p;
    cout << primitive_root(p) << "\n";
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