// @brief Bernoulli Number
#define PROBLEM "https://judge.yosupo.jp/problem/bernoulli_number"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/math/poly/series.hpp"
#include "cp-algo/math/poly/transform.hpp"
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

// EGF of Bk
polyn bernoulli(size_t n) {
    return inv((expx<base>(n+1) - polyn(1)).div_xk(1), n);
}

void solve() {
    int n;
    cin >> n;
    invborel(bernoulli(n+1)).print(n+1);
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
