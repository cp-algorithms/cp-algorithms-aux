// @brief Ordered Set
#define PROBLEM "https://judge.yosupo.jp/problem/ordered_set"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <bits/stdc++.h>
//#include "blazingio/blazingio.min.hpp"
#include "cp-algo/structures/fenwick_set.hpp"
#include "cp-algo/util/compress_coords.hpp"

using namespace std;
using cp_algo::structures::fenwick_set;

void solve() {
    int n, q;
    cin >> n >> q;
    cp_algo::big_vector<int> a(n, 0);
    cp_algo::big_vector<reference_wrapper<int>> coords;
    for(auto &it: a) {
        cin >> it;
        coords.push_back(ref(it));
    }
    cp_algo::big_vector<pair<int, int>> queries(q, pair{0, 0});
    for(auto &[t, x]: queries) {
        cin >> t >> x;
        if(t != 2) {
            coords.push_back(ref(x));
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
            auto res = (int)me.find_by_order(x-1);
            cout << (res == -1 ? -1 : values[res]) << '\n';
        } else if(t == 3) {
            cout << me.order_of_key(x+1) << '\n';
        } else if(t == 4) {
            auto res = (int)me.pre_upper_bound(x);
            cout << (res == -1 ? -1 : values[res]) << '\n';
        } else if(t == 5) {
            auto res = (int)me.lower_bound(x);
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
