// @brief Convolution mod $10^9+7$
#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod_1000000007"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("tune=native")
#include "cp-algo/math/fft.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 1e9 + 7;
using base = modint<mod>;

void solve() {
    int n, m;
    cin >> n >> m;
    vector<base> a(n), b(m);
    copy_n(istream_iterator<base>(cin), n, begin(a));
    copy_n(istream_iterator<base>(cin), m, begin(b));
    fft::mul(a, b);
    ranges::copy(a, ostream_iterator<base>(cout, " "));
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
