// @brief Convolution GF(2^64)
#define PROBLEM "https://judge.yosupo.jp/problem/convolution_F_2_64"
#pragma GCC optimize("O3,unroll-loops")
#include <bits/allocator.h>
#pragma GCC target("avx2,vpclmulqdq")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#define CP_ALGO_CHECKPOINT
#include "cp-algo/util/big_alloc.hpp"
#include "cp-algo/number_theory/nimber.hpp"
#include "cp-algo/math/karatsuba.hpp"
#include <bits/stdc++.h>

using namespace std;

void solve() {
    size_t n, m;
    cin >> n >> m;
    cp_algo::big_vector<f2_64> a(n), b(m);
    for (auto &x : a) cin >> x.r;
    for (auto &x : b) cin >> x.r;
    auto c = cp_algo::math::karatsuba(a, b);
    for (auto x : c) cout << x.r << " ";
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
