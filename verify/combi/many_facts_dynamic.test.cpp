// @brief Many Factorials (dynamic Mod)
#define PROBLEM "https://judge.yosupo.jp/problem/many_factorials"
#pragma GCC optimize("O3,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/factorials.hpp"

using namespace std;
using base = cp_algo::math::dynamic_modint<>;

void solve() {
    int n;
    cin >> n;
    base::switch_mod(998'244'353);
    vector<base> args(n);
    for(auto &x : args) {cin >> x;}
    cp_algo::checkpoint("read");
    auto res = facts<true, 100'000>(args);
    for(auto it: res) {cout << it << "\n";}
    cp_algo::checkpoint("write");
    cp_algo::checkpoint<1>();
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
