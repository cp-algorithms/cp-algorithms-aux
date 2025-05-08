// @brief Longest Increasing Subsequence
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/longest_increasing_subsequence
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/structures/fenwick.hpp"
#include "cp-algo/util/compress_coords.hpp"
#include <bits/stdc++.h>

using namespace std;

void solve() {
    int n;
    cin >> n;
    vector<reference_wrapper<int>> coords;
    vector<int> a(n);
    for(auto &it: a) {
        cin >> it;
        coords.push_back(ref(it));
    }
    auto b = cp_algo::compress_coords(coords);
    struct longest {
        int len = 0, last = 0;
        auto operator <=>(longest const&) const = default;
    };
    cp_algo::structures::fenwick_max me(vector<longest>(n + 1));
    vector<int> pre(n);
    for(auto [i, it]: a | views::enumerate) {
        auto [len, last] = me.prefix_fold(it);
        me.update(it, {len + 1, (int)i});
        pre[i] = last;
    }
    auto [len, last] = me.prefix_fold(n);
    cout << len << '\n';
    vector<int> ans(len);
    while(len--) {
        ans[len] = last;
        last = pre[last];
    }
    for(auto it: ans) {
        cout << it << ' ';
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
 