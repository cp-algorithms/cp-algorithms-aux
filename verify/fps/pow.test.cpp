// @brief Lazy FPS: Pow of Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/pow_of_formal_power_series"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/fps.hpp"

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n;
    int64_t m;
    cin >> n >> m;
    polyn::Vector a(n);
    for(auto &it: a) {cin >> it;}
    pow(fps<base>(polyn(std::move(a))), m).prefix(n).print(n);
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
