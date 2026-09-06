// @brief Composition of Formal Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/composition_of_formal_power_series"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/poly/compose.hpp"

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n;
    int64_t m;
    cin >> n;
    polyn::Vector a(n), b(n);
    for(auto &it: a) {cin >> it;}
    for(auto &it: b) {cin >> it;}
    compose(polyn(std::move(a)), polyn(std::move(b)), n).print(n);
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
