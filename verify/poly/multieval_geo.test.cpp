// @brief Multipoint Evaluation (Geometric Sequence)
#define PROBLEM "https://judge.yosupo.jp/problem/multipoint_evaluation_on_geometric_sequence"
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
    int n, m, a, r;
    cin >> n >> m >> a >> r;
    polyn::Vector f(n);
    for(auto &it: f) {cin >> it;}
    chirpz(mulx(polyn(std::move(f)), a), r, m).print(m);
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
