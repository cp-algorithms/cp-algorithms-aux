// @brief Lazy FPS: Log of Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/log_of_formal_power_series"
#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/fps.hpp"

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n;
    cin >> n;
    polyn::Vector a(n);
    for(auto &x: a) {cin >> x;}
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
