// @brief Polynomial Interpolation
#define PROBLEM "https://judge.yosupo.jp/problem/polynomial_interpolation"
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/io.hpp"
#include "cp-algo/math/poly/eval.hpp"

using cp_algo::io::read_reduced;

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n;
    cin >> n;
    polyn::Vector x(n), y(n);
    read_reduced(x);
    read_reduced(y);
    inter(x, y).print(n);
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
