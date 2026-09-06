// @brief Consecutive Terms of Linear Recursion
#define PROBLEM "https://judge.yosupo.jp/problem/consecutive_terms_of_linear_recurrent_sequence"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/poly/recurrence.hpp"

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int d, M;
    int64_t k;
    cin >> d >> k >> M;
    polyn::Vector a(d), c(d);
    for(auto &it: a) {cin >> it;}
    for(auto &it: c) {cin >> it;}
    polyn A = polyn(std::move(a));
    polyn Q = polyn::xk(0) - polyn(c).mul_xk(1);
    polyn P = (A * Q).mod_xk(d);
    (P * inv(Q, k - d, M + d)).div_xk(d).print(M);
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