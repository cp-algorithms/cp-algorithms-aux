// @brief Dynamic Range Affine Range Sum
#define PROBLEM "https://judge.yosupo.jp/problem/dynamic_sequence_range_affine_range_sum"
#include "cp-algo/algebra/modular.hpp"
#include "cp-algo/data_structures/treap/metas/reverse.hpp"
#include "cp-algo/data_structures/treap.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::data_structures;

using base = cp_algo::algebra::modular<998244353>;
using meta = treap::metas::reverse_meta<base>;
using node = treap_node<meta>;
using treap_t = node::treap;

void solve() {
    istream_iterator<int> input(cin);
    int n = *input++;
    int q = *input++;
    vector<treap_t> nodes(n);
    generate_n(begin(nodes), n, [&](){
        return new node(meta(*input++));
    });
    auto me = node::build(nodes);

    while(q--) {
        int t = *input++;
        if(t == 0) {
            int i = *input++;
            base x = *input++;
            node::insert(me, i, new node(meta(x)));
        } else if(t == 1) {
            node::erase(me, *input++);
        } else if(t == 2) {
            int l = *input++;
            int r = *input++;
            node::exec_on_segment(me, l, r, [](auto &t) {
                _safe_meta(t, reverse = 1);
            });
        } else if(t == 3) {
            int l = *input++;
            int r = *input++;
            base b = *input++;
            base c = *input++;
            node::exec_on_segment(me, l, r, [b, c](auto &t) {
                _safe_meta(t, add_push(meta::lin(b, c)));
            });
        } else {
            int l = *input++;
            int r = *input++;
            node::exec_on_segment(me, l, r, [](auto t) {
                cout << _safe_meta(t, sum) << "\n";
            });
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
