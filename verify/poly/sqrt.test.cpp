// @brief Sqrt of Power Series
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/sqrt_of_formal_power_series
#include "cp-algo/math/poly.hpp"
#include "cp-algo/number_theory/modint.hpp"
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
    auto res = polyn(a).sqrt(n);
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
