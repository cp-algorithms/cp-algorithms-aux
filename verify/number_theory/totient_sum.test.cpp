// @brief Sum of Totient Function
#define PROBLEM "https://judge.yosupo.jp/problem/sum_of_totient_function"
#pragma GCC optimize("Ofast,unroll-loops")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/number_theory/modint.hpp"
#include "cp-algo/number_theory/dirichlet.hpp"
#include "cp-algo/util/big_alloc.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;
using base = modint<998244353>;

void solve() {
    int64_t n;
    cin >> n;
    auto H = floors(n) | views::transform([](int64_t k) {
        return base(k) * base(k + 1) / 2;
    }) | ranges::to<vector>();
    auto G = floors(n) | views::transform([](int64_t k) {
        return base(k);
    }) | ranges::to<vector>();
    auto F = Dirichlet_div(H, G, n);
    cout << F.back() << "\n";
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    //cin >> t;
    while(t--) {
        solve();
    }
}
