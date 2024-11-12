// @brief Ordered Set
#define PROBLEM "https://judge.yosupo.jp/problem/ordered_set"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/structures/fenwick_set.hpp"
#include "cp-algo/util/compress_coords.hpp"
#include <bits/stdc++.h>

using namespace std;
using cp_algo::structures::fenwick_set;

void solve() {
    int n, q;
    cin >> n >> q;
    vector a(n, 0);
    vector<int*> coords;
    for(auto &it: a) {
        cin >> it;
        coords.push_back(&it);
    }
    vector queries(q, pair{0, 0});
    for(auto &[t, x]: queries) {
        cin >> t >> x;
        if(t != 2) {
            coords.push_back(&x);
        }
    }
    auto values = cp_algo::compress_coords(coords);
    const int maxc = 1e6;
    fenwick_set<maxc> me(a);
    for(auto [t, x]: queries) {
        if(t == 0) {
            me.insert(x);
        } else if(t == 1) {
            me.erase(x);
        } else if(t == 2) {
            int res = me.find_by_order(x-1);
            cout << (res == -1 ? -1 : values[res]) << '\n';
        } else if(t == 3) {
            cout << me.order_of_key(x+1) << '\n';
        } else if(t == 4) {
            int res = me.pre_upper_bound(x);
            cout << (res == -1 ? -1 : values[res]) << '\n';
        } else if(t == 5) {
            int res = me.lower_bound(x);
            cout << (res == -1 ? -1 : values[res]) << '\n';
        }
    }
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
