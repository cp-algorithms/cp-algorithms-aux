// @brief Polynomial Interpolation (Geometric Sequence)
#define PROBLEM "https://judge.yosupo.jp/problem/polynomial_interpolation_on_geometric_sequence"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/poly/chirpz.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n, a, r;
    cin >> n >> a >> r;
    polyn::Vector y(n);
    for(auto &it: y) {cin >> it;}
    mulx(chirpz_inverse(polyn(std::move(y)), r, n), base(a).inv()).print(n);
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    while(t--) {
        solve();
    }
}
