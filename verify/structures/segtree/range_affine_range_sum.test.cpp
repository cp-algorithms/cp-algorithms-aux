// @brief Range Affine Range Sum
#define PROBLEM "https://judge.yosupo.jp/problem/range_affine_range_sum"
#include "cp-algo/structures/segtree/metas/affine.hpp"
#include "cp-algo/structures/segtree.hpp"
#include "cp-algo/number_theory/modint.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::structures;

using base = cp_algo::math::modint<998244353>;
using meta = segtree::metas::affine_meta<base>;

void solve() {
    int n, q;
    cin >> n >> q;
    cp_algo::big_vector<meta> a(n);
    for(int i = 0; i < n; i++) {
        int ai;
        cin >> ai;
        a[i] = {ai};
    }
    segtree_t<meta> me(a);
    while(q--) {
        int t;
        cin >> t;
        if(t == 0) {
            int l, r, b, c;
            cin >> l >> r >> b >> c;
            me.exec_on_segment(l, r, [&](auto& meta) {
                meta.to_push.prepend(meta::lin(b, c));
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
