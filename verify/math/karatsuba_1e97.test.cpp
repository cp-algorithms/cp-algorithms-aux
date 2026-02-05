// @brief Convolution (Mod 1,000,000,007, Karatsuba)
#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod_1000000007"
#pragma GCC optimize("O3,unroll-loops")
#include <bits/allocator.h>
#pragma GCC target("avx2,vpclmulqdq")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#define CP_ALGO_CHECKPOINT
#include "cp-algo/math/karatsuba.hpp"
#include <bits/stdc++.h>

using namespace std;

const int mod = 1e9+7;
using base = cp_algo::math::modint<mod>;

void solve() {
    size_t n, m;
    cin >> n >> m;
    cp_algo::big_vector<base> a(n), b(m);
    for (auto &x : a) cin >> x;
    for (auto &x : b) cin >> x;
    auto c = cp_algo::math::karatsuba(a, b);
    for (auto x : c) cout << x << " ";
    cout << "\n";
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    while(t--) {
        solve();
    }
}
