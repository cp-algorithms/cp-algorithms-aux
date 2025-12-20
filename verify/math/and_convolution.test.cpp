// @brief Bitwise And Convolution
#define PROBLEM "https://judge.yosupo.jp/problem/bitwise_and_convolution"
#pragma GCC optimize("O3,unroll-loops")
#include <bits/allocator.h>
#pragma GCC target("avx2")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#define CP_ALGO_CHECKPOINT
#include "cp-algo/number_theory/modint.hpp"
#include "cp-algo/util/big_alloc.hpp"
#include "cp-algo/util/checkpoint.hpp"
#include "cp-algo/math/and_convolution.hpp"
#include <bits/stdc++.h>

using namespace std;

const int mod = 998244353;
using base = cp_algo::math::modint<mod>;

void solve() {
    uint32_t n;
    cin >> n;
    uint32_t N = 1u << n;
    cp_algo::big_vector<base> a(N), b(N);
    for (auto &it : a) {cin >> it;}
    for (auto &it : b) {cin >> it;}
    cp_algo::checkpoint("read");
    cp_algo::math::and_convolution_inplace(a, b);
    for (auto it : a) {cout << it << ' ';}
    cp_algo::checkpoint("write");
    cp_algo::checkpoint<1>();
}

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    while (t--) {
        solve();
    }
}
