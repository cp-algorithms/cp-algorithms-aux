// @brief Range Reverse Range Sum
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/range_reverse_range_sum
#include "cp-algo/structures/treap/metas/reverse.hpp"
#include "cp-algo/structures/treap.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::structures::treap;

using meta = metas::reverse_meta<int64_t>;
using node_t = node<meta>;
using treap = node_t::treap;

void solve() {
    istream_iterator<int> input(cin);
    int n = *input++;
    int q = *input++;
    vector<treap> nodes(n);
    generate_n(begin(nodes), n, [&](){
        return node_t::make_treap(meta(*input++));
    });
    auto me = node_t::build(nodes);
    while(q--) {
        int t = *input++;
        int l = *input++;
        int r = *input++;
        if(t == 0) {
            node_t::exec_on_segment(me, l, r, [](auto &t) {
                _safe_meta(t, reverse = true);
            });
        } else {
            node_t::exec_on_segment(me, l, r, [](auto const& t) {
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
