// @brief Dirichlet Convolution and Prefix Sums
#define PROBLEM "https://judge.yosupo.jp/problem/dirichlet_convolution_and_prefix_sums"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <iostream>
//#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/big_alloc.hpp"
#include "cp-algo/number_theory/modint.hpp"
#include "cp-algo/number_theory/dirichlet.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;
using base = modint<998244353>;

void solve() {
    int64_t n;
    cin >> n;
    auto [_, m] = floor_stats(n);
    cp_algo::big_vector<base> F(m+1), G(m+1);
    for (int i = 1; i <= m; ++i) {
        cin >> F[i];
    }
    for (int i = 1; i <= m; ++i) {
        cin >> G[i];
    }
    auto H = Dirichlet_mul(F, G, n);
    for (int i = 1; i <= m; ++i) {
        cout << H[i] << " \n"[i == m];
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
