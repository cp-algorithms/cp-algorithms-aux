// @brief Point Add Range Sum
#define PROBLEM "https://judge.yosupo.jp/problem/point_add_range_sum"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/structures/fenwick.hpp"
#include <bits/stdc++.h>

using namespace std;

void solve() {
    int n, q;
    cin >> n >> q;
    vector<int64_t> a(n + 1);
    for(auto &it: a | views::drop(1)) {cin >> it;}
    cp_algo::structures::fenwick<int64_t> me(move(a));
    for(int i = 0; i < q; i++) {
        int t, x, y;
        cin >> t >> x >> y;
        if(t == 0) {
            me.add(x, y);
        } else {
            cout << me.range_sum(x, y) << '\n';
        }
    }
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
