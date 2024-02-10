// @brief Range Affine Range Sum
#define PROBLEM "https://judge.yosupo.jp/problem/range_affine_range_sum"
#include "cp-algo/data_structures/metas/affine.hpp"
#include "cp-algo/data_structures/segment_tree.hpp"
#include "cp-algo/algebra/modular.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace data_structures::segment_tree;

using base = algebra::modular<998244353>;
using meta = metas::affine_meta<base>;

void solve() {
    int n, q;
    cin >> n >> q;
    vector<meta> a(n);
    for(int i = 0; i < n; i++) {
        int ai;
        cin >> ai;
        a[i] = {ai};
    }
    segment_tree<meta> me(a);
    while(q--) {
        int t;
        cin >> t;
        if(t == 0) {
            int l, r, b, c;
            cin >> l >> r >> b >> c;
            me.exec_on_segment(l, r, [&](auto& meta) {
                meta.to_push = meta::lin(b, c) * meta.to_push;
            });
        } else {
            int l, r;
            cin >> l >> r;
            base ans = 0;
            me.exec_on_segment(l, r, [&](auto meta) {
                ans += meta.sum;
            });
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
