// @brief Inv of Polynomials
#define PROBLEM "https://judge.yosupo.jp/problem/inv_of_polynomials"
#include "cp-algo/math/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n, m;
    cin >> n >> m;
    vector<base> a(n), b(m);
    copy_n(istream_iterator<base>(cin), n, begin(a));
    copy_n(istream_iterator<base>(cin), m, begin(b));
    auto res = polyn(a).inv_mod(polyn(b));
    if(res) {
        cout << res->deg() + 1 << endl;
        res->print();
    } else {
        cout << -1 << endl;
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