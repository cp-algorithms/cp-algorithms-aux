// @brief Dirichlet Inverse and Prefix Sums
#define PROBLEM "https://judge.yosupo.jp/problem/dirichlet_inverse_and_prefix_sums"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/number_theory/modint.hpp"
#include "cp-algo/number_theory/dirichlet.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;
using base = modint<998244353>;

void solve() {
    int64_t n;
    cin >> n;
    auto [_, m] = floor_stats(n);
    vector<base> H(m+1), G(m+1);
    for (int i = 1; i <= m; ++i) {
        cin >> G[i];
    }
    for (int i = 1; i <= m; ++i) {
        H[i] = 1;
    }
    auto F = Dirichlet_div(H, G, n);
    for (int i = 1; i <= m; ++i) {
        cout << F[i] << " \n"[i == m];
    }
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    cin >> t;
    while(t--) {
        solve();
    }
}
