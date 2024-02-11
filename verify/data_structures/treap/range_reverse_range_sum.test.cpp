// @brief Range Reverse Range Sum
#define PROBLEM "https://judge.yosupo.jp/problem/range_reverse_range_sum"
#include "cp-algo/data_structures/treap/metas/reverse.hpp"
#include "cp-algo/data_structures/treap.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::data_structures;

using meta = treap::metas::reverse_meta<int64_t>;
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
        int l = *input++;
        int r = *input++;
        if(t == 0) {
            node::exec_on_segment(me, l, r, [](auto &t) {
                _safe_meta(t, reverse = true);
            });
        } else {
            node::exec_on_segment(me, l, r, [](auto const& t) {
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
