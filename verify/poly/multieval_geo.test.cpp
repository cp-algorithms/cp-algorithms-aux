// @brief Multipoint Evaluation (Geometric Sequence)
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/multipoint_evaluation_on_geometric_sequence
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/math/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n, m, a, r;
    cin >> n >> m >> a >> r;
    polyn::Vector f(n);
    copy_n(istream_iterator<base>(cin), n, begin(f));
    polyn(f).mulx(a).chirpz(r, m).print(m);
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
