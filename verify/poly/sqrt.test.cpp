// @brief Sqrt of Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/sqrt_of_formal_power_series"
#pragma GCC optimize("O3,unroll-loops")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/io.hpp"
#include "cp-algo/math/poly/series/sqrt.hpp"
#include "cp-algo/math/poly/series.hpp"

using cp_algo::io::read_reduced;

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n;
    cin >> n;
    polyn::Vector a(n);
    read_reduced(a);
    auto res = sqrt(polyn(std::move(a)), n);
    if(res) {
        res->print(n);
    } else {
        cout << -1 << "\n";
    }
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
