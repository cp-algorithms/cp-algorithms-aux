// @brief Polynomial Composite Set Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/polynomial_composite_set_power_series"
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
    size_t M, n;
    cin >> M >> n;
    size_t N = 1 << n;
    cp_algo::big_vector<base> f(M);
    for(auto &it: f) {cin >> it;}
    cp_algo::big_vector<base> a(N);
    for(auto &it: a) {cin >> it;}
    cp_algo::checkpoint("read");
    auto c = cp_algo::math::subset_compose<base>(f, a);
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
