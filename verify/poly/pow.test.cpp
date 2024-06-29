// @brief Pow of Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/pow_of_formal_power_series"
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
    int64_t m;
    cin >> n >> m;
    vector<base> a(n);
    copy_n(istream_iterator<base>(cin), n, begin(a));
    polyn(a).pow(m, n).print(n);
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
