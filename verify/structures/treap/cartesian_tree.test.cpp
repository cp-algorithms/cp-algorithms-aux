// @brief Build Cartesian Tree
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/cartesian_tree
#include "cp-algo/structures/treap/metas/base.hpp"
#include "cp-algo/structures/treap.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::structures::treap;

struct val_meta: metas::base_meta {
    int val;
    val_meta(int val): val(val){}
};

using node_t = node<val_meta>;
using treap = node_t::treap;

void solve() {
    istream_iterator<int> input(cin);
    int n = *input++;
    vector<treap> nodes(n);
    for(int i = 0; i < n; i++) {
        nodes[i] = node_t::make_treap(val_meta(i), *input++);
    }
    auto me = node_t::build(nodes);
    vector<int> p(n, -1);
    node_t::exec_on_each(me, [&](auto t) {
        for(auto child: t->children) {
            if(child) {
                p[_safe_meta(child, val)] = _safe_meta(t, val);
            }
        }
    });
    for(int i = 0; i < n; i++) {
        cout << (p[i] == -1 ? i : p[i]) << ' ';
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
