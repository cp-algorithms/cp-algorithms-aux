// @brief Characteristic Polynomial
#define PROBLEM "https://judge.yosupo.jp/problem/characteristic_polynomial"
#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_MAXN 1 << 10
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
    cin >> n;
    matrix<base> A(n);
    A.read();
    auto blocks = frobenius_form(A);
    reduce(begin(blocks), end(blocks), polyn(1), multiplies{}).print();
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
