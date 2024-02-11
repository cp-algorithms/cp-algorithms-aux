// @brief Inv of Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/inv_of_formal_power_series"
#include "cp-algo/algebra/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::algebra;

const int mod = 998244353;
using base = modular<mod>;
using polyn = poly_t<base>;

void solve() {
    int n;
    cin >> n;
    vector<base> a(n);
    copy_n(istream_iterator<base>(cin), n, begin(a));
    polyn(a).inv(n).print(n);
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
