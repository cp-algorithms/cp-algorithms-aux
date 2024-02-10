#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod_1000000007"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("avx2,tune=native")
#include "cp-algo/algebra/polynomial.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace algebra;

const int mod = 1e9 + 7;
typedef modular<mod> base;
typedef poly<base> polyn;

void solve() {
    int n, m;
    cin >> n >> m;
    vector<base> a(n), b(m);
    copy_n(istream_iterator<base>(cin), n, begin(a));
    copy_n(istream_iterator<base>(cin), m, begin(b));
    (polyn(a) * polyn(b)).print(n + m - 1);
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t;
    t = 1;// cin >> t;
    while(t--) {
        solve();
    }
}