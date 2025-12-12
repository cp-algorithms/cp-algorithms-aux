// @brief Factorize
#define PROBLEM "https://judge.yosupo.jp/problem/factorize"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/number_theory/factorize.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

void solve() {
    int64_t m;
    cin >> m;
    auto res = to<vector>(factorize(m));
    ranges::sort(res);
    cout << size(res) << " ";
    ranges::copy(res, ostream_iterator<int64_t>(cout, " "));
    cout << "\n";
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
