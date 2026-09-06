// @brief Lazy FPS: Log of Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/log_of_formal_power_series"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/math/fps.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n;
    cin >> n;
    polyn::Vector a(n);
    copy_n(istream_iterator<base>(cin), n, begin(a));
    log(fps<base>(polyn(std::move(a)))).prefix(n).print(n);
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
