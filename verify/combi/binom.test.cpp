// @brief Binomial Coefficient (Prime Mod)
#define PROBLEM "https://judge.yosupo.jp/problem/binomial_coefficient_prime_mod"
#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_MAXN 1e7
#include "cp-algo/number_theory/modint.hpp"
#include "cp-algo/math/combinatorics.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo;
using namespace math;
using base = dynamic_modint<>;

void solve() {
    int n, r;
    cin >> n >> r;
    cout << binom<base>(n, r) << "\n";
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    cin >> t;
    int m;
    cin >> m;
    base::switch_mod(m);
    while(t--) {
        solve();
    }
}