// @brief Polynomial Taylor Shift
#define PROBLEM "https://judge.yosupo.jp/problem/polynomial_taylor_shift"
#pragma GCC optimize("O3,unroll-loops")
#include "cp-algo/math/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n, c;
    cin >> n >> c;
    polyn::Vector a(n);
    copy_n(istream_iterator<base>(cin), n, begin(a));
    polyn(a).shift(c).print(n);
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
