// @brief Multipoint Evaluation
#define PROBLEM "https://judge.yosupo.jp/problem/multipoint_evaluation"
#pragma GCC optimize("O3,unroll-loops")
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
    polyn::Vector f(n), x(m);
    copy_n(istream_iterator<base>(cin), n, begin(f));
    copy_n(istream_iterator<base>(cin), m, begin(x));
    polyn(polyn(f).eval(x)).print(m);
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
