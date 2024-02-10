// @brief Sqrt of Formal Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/sqrt_of_formal_power_series"
#include "cp-algo/algebra/polynomial.hpp"
#include "cp-algo/algebra/modular.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::algebra;

const int mod = 998244353;
typedef modular<mod> base;
typedef poly<base> polyn;

void solve() {
    int n;
    cin >> n;
    vector<base> a(n);
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
