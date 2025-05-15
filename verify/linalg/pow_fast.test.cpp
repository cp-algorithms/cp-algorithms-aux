// @brief Pow of Matrix (Frobenius)
#define PROBLEM "https://judge.yosupo.jp/problem/pow_of_matrix"
#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_MAXN 256
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/linalg/frobenius.hpp"

using namespace std;
using namespace cp_algo::math;
using namespace cp_algo::linalg;

const int64_t mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    size_t n;
    uint64_t k;
    cin >> n >> k;
    matrix<base> A(n);
    A.read();
    frobenius_pow(A, k).print();
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    //cin >> t;
    while(t--) {
        solve();
    }
}
