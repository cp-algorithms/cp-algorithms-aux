// @brief Point Add Range Sum
#define PROBLEM "https://judge.yosupo.jp/problem/point_add_range_sum"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/structures/fenwick.hpp"
#include <bits/stdc++.h>

using namespace std;

void solve() {
    int n, q;
    cin >> n >> q;
    cp_algo::big_vector<int64_t> a(n + 1);
    for(auto &it: a | views::drop(1)) {cin >> it;}
    cp_algo::structures::fenwick me(move(a));
    for(int i = 0; i < q; i++) {
        int t, x, y;
        cin >> t >> x >> y;
        if(t == 0) {
            me.update(x, y);
        } else {
            cout << me.range_fold(x, y) << '\n';
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
