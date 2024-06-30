// @brief Kth term of Linearly Recurrent Sequence
#define PROBLEM "https://judge.yosupo.jp/problem/kth_term_of_linearly_recurrent_sequence"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/math/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int64_t d, k;
    cin >> d >> k;
    vector<base> a(d), c(d);
    copy_n(istream_iterator<base>(cin), d, begin(a));
    copy_n(istream_iterator<base>(cin), d, begin(c));
    polyn Q = polyn(1) - polyn(c).mul_xk(1);
    polyn P = (polyn(a) * Q).mod_xk(d);
    cout << polyn::kth_rec(P, Q, k) << endl;
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
