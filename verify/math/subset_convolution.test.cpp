// @brief Subset Convolution
#define PROBLEM "https://judge.yosupo.jp/problem/subset_convolution"
#pragma GCC optimize("O3,unroll-loops")
#include <bits/allocator.h>
#pragma GCC target("avx2")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#define CP_ALGO_CHECKPOINT
#include "cp-algo/number_theory/modint.hpp"
#include "cp-algo/math/subset_convolution.hpp"
#include <bits/stdc++.h>

using namespace std;

const int mod = 998244353;
using base = cp_algo::math::modint<mod>;

void solve() {
    uint32_t n;
    cin >> n;
    size_t N = 1 << n;
    cp_algo::big_vector<base> a(N), b(N);
    for(auto &it: a) {cin >> it;}
    for(auto &it: b) {cin >> it;}
    cp_algo::checkpoint("read");
    auto c = cp_algo::math::subset_convolution<base>(a, b);
    for(auto &it: c) {cout << it << ' ';}
    cp_algo::checkpoint("write");
    cp_algo::checkpoint<1>();
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t;
    t = 1;// cin >> t;
    while(t--) {
        solve();
    }
}
