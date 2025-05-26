// @brief Many Factorials
#define PROBLEM "https://judge.yosupo.jp/problem/many_factorials"
#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/factorials.hpp"

using namespace std;
using base = cp_algo::math::modint<998244353>;

void solve() {
    int n;
    cin >> n;
    vector<base> args(n);
    for(auto &x : args) {cin >> x;}
    cp_algo::checkpoint("read");
    auto res = facts(args);
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
