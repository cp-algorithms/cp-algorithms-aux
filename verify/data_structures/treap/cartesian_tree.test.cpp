// @brief Build Cartesian Tree
#define PROBLEM "https://judge.yosupo.jp/problem/cartesian_tree"
#include "cp-algo/data_structures/treap/metas/base.hpp"
#include "cp-algo/data_structures/treap.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::data_structures;

struct val_meta: treap::metas::base_meta {
    int val;
    val_meta(int val): val(val){}
};

using node = treap_node<val_meta>;
using treap_t = node::treap;

void solve() {
    istream_iterator<int> input(cin);
    int n = *input++;
    vector<treap_t> nodes(n);
    for(int i = 0; i < n; i++) {
        nodes[i] = new node(val_meta(i), *input++);
    }
    auto me = node::build(nodes);
    vector<int> p(n, -1);
    node::exec_on_each(me, [&](auto t) {
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
