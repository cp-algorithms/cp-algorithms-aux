// @brief Bell Number
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/bell_number
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/math/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n;
    cin >> n;
    (polyn::expx(n+1) - polyn(1)).exp_inplace(n+1).invborel().print(n+1);
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
