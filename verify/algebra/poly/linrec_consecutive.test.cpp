// @brief Consecutive Terms of Linear Recursion
#define PROBLEM "https://judge.yosupo.jp/problem/consecutive_terms_of_linear_recurrent_sequence"
#include "cp-algo/algebra/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::algebra;

const int mod = 998244353;
using base = modular<mod>;
using polyn = poly_t<base>;

void solve() {
    int d, M;
    int64_t k;
    cin >> d >> k >> M;
    vector<base> a(d), c(d);
    copy_n(istream_iterator<base>(cin), d, begin(a));
    copy_n(istream_iterator<base>(cin), d, begin(c));
    polyn A = polyn(a);
    polyn Q = polyn::xk(0) - polyn(c).mul_xk(1);
    polyn P = (A * Q).mod_xk(d);
    A = (P * Q.inv(2 * d)).mod_xk(2 * d);
    A = A.reverse(2 * d);
    polyn shift = polyn::xk(1).powmod(k, Q.reverse(d + 1));
    polyn An = (A * shift).div_xk(d - 1).mod_xk(d + 1).reverse(d + 1);
    polyn Pn = (An * Q).mod_xk(d);
    (Pn * Q.inv(M)).print(M);
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