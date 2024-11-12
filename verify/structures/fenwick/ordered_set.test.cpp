// @brief Ordered Set
#define PROBLEM "https://judge.yosupo.jp/problem/ordered_set"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/structures/fenwick_set.hpp"
#include <bits/stdc++.h>

using namespace std;
using cp_algo::structures::fenwick_set;

void solve() {
    int n, q;
    cin >> n >> q;
    vector a(n, 0);
    for(auto &it: a) {cin >> it;}
    auto nums = a;
    vector queries(q, pair{0, 0});
    for(auto &[t, x]: queries) {
        cin >> t >> x;
        if(t == 0) {
            nums.push_back(x);
        }
    }
    nums.push_back(0);
    ranges::sort(nums);
    nums.erase(ranges::unique(nums).begin(), end(nums));
    auto hashify = [&](int x) {
        return ranges::lower_bound(nums, x) - begin(nums);
    };
    ranges::transform(a, begin(a), hashify);
    const int maxc = 1e6;
    fenwick_set<maxc> me(a);
    for(auto [t, x]: queries) {
        if(t == 0) {
            me.insert(hashify(x));
        } else if(t == 1) {
            int y = hashify(x);
            if(nums[y] == x) {
                me.erase(y);
            }
        } else if(t == 2) {
            int res = me.find_by_order(x-1);
            cout << (res == -1 ? -1 : nums[res]) << '\n';
        } else if(t == 3) {
            cout << me.order_of_key(hashify(x+1)) << '\n';
        } else if(t == 4) {
            int res = me.pre_upper_bound(hashify(x+1)-1);
            cout << (res == -1 ? -1 : nums[res]) << '\n';
        } else if(t == 5) {
            int res = me.lower_bound(hashify(x));
            cout << (res == -1 ? -1 : nums[res]) << '\n';
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
