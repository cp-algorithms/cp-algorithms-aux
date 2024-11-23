// @brief Shift of Sampling Points of Polynomial
#define PROBLEM "https://judge.yosupo.jp/problem/shift_of_sampling_points_of_polynomial"
#define CP_ALGO_MAXN 1 << 20
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/math/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

// TODO: Use single-convolution approach
void solve() {
    int n, m, c;
    cin >> n >> m >> c;
    vector<base> a(n);
    copy_n(istream_iterator<base>(cin), n, begin(a));
    polyn A = polyn(a);
    polyn Q = polyn({1, -1}).pow(n, n + 1);
    A -= ((A * Q).div_xk(n).mod_xk(m) * Q.inv(m)).mod_xk(m).mul_xk(n);
    A = A.reverse(n + m);
    polyn shift = polyn({1, -1}).pow(c, n).shift(1).mulx(-1);
    auto R = (A.div_xk(n - 1) * shift) + (A.mod_xk(n - 1) * shift).div_xk(n - 1);
    R.div_xk(1).reverse(m).print(m);
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
