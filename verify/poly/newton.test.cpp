// @brief Conversion to Newton Basis
#define PROBLEM "https://judge.yosupo.jp/problem/conversion_from_monomial_basis_to_newton_basis"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/math/poly/eval.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n;
    cin >> n;
    polyn::Vector a(n), p(n);
    copy_n(istream_iterator<base>(cin), n, begin(a));
    copy_n(istream_iterator<base>(cin), n, begin(p));
    to_newton(polyn(std::move(a)), p).print(n);
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
