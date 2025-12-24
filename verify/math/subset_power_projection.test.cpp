// @brief Power Projection of Set Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/power_projection_of_set_power_series"
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
    uint32_t n, M;
    cin >> n >> M;
    size_t N = 1 << n;
    cp_algo::big_vector<base> g(N);
    for(auto &it: g) {cin >> it;}
    cp_algo::big_vector<base> w(N);
    for(auto &it: w) {cin >> it;}
    cp_algo::checkpoint("read");
    auto c = cp_algo::math::subset_power_projection<base>(g, w, M);
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
