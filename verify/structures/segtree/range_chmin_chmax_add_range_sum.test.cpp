// @brief Range Chmin Chmax Add Range Sum
#define PROBLEM "https://judge.yosupo.jp/problem/range_chmin_chmax_add_range_sum"
#include "cp-algo/structures/segtree/metas/chmin_chmax_add.hpp"
#include "cp-algo/structures/segtree.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::structures;
using meta = segtree::metas::chmin_chmax_sum_meta;

void solve() {
    int n, q;
    cin >> n >> q;
    vector<meta> a(n);
    for(int i = 0; i < n; i++) {
        int64_t ai;
        cin >> ai;
        a[i] = {ai};
    }
    segtree_t<meta> me(a);
    while(q--) {
        int t, l, r;
        int64_t b;
        cin >> t >> l >> r;
        if(t == 0) {
            cin >> b;
            me.exec_on_segment(l, r,
                [b](auto& meta) {meta.chmin = b;},
                meta::proceed_chmin(b), meta::stop_chmin(b));
        } else if(t == 1) {
            cin >> b;
            me.exec_on_segment(l, r,
                [b](auto& meta) {meta.chmax = b;},
                meta::proceed_chmax(b), meta::stop_chmax(b));
        } else if(t == 2) {
            cin >> b;
            me.exec_on_segment(l, r,
                [b](auto& meta) {meta.add = b;});
        } else {
            int64_t ans = 0;
            me.exec_on_segment(l, r, [&](auto& meta) {
                ans += meta.sum;});
            cout << ans << "\n";
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
