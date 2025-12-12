// @brief Enumerate Quotients
#define PROBLEM "https://judge.yosupo.jp/problem/enumerate_quotients"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/number_theory/dirichlet.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

void solve() {
    int64_t n;
    cin >> n;
    auto res = floors(n) | views::drop(1);
    cout << size(res) << "\n";
    for(auto it: res) {
        cout << it << " ";
    }
    cout << "\n";
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    //cin >> t;
    while(t--) {
        solve();
    }
}
