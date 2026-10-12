#include <bits/stdc++.h>
#ifdef TEST_MINIFIED
#include "cp-algo/min/structures/segtree/metas/chmin_chmax_add.hpp"
#include "cp-algo/min/structures/segtree.hpp"
#else
#include "cp-algo/structures/segtree/metas/chmin_chmax_add.hpp"
#include "cp-algo/structures/segtree.hpp"
#endif
using meta = cp_algo::structures::segtree::metas::chmin_chmax_sum_meta;
int main() {
    std::mt19937 rng(7231);
    for(int test = 0; test < 200; ++test) {
        int n = 1 + rng() % 30;
        std::vector<int64_t> expected(n);
        cp_algo::structures::segtree_t<meta> tree(cp_algo::big_vector<meta>(n, meta(0)));
        for(int q = 0; q < 300; ++q) {
            int l = rng() % n, r = l + 1 + rng() % (n - l);
            int64_t b = q % 3 == 0 ? -meta::inf : q % 3 == 1 ? meta::inf : int64_t(rng() % 101) - 50;
            if(q % 2 == 0) {
                tree.exec_on_segment(l, r, [b](auto& m) {m.chmin = b;},
                                     meta::proceed_chmin(b), meta::stop_chmin(b));
                for(int i = l; i < r; ++i) expected[i] = std::min(expected[i], b);
            } else {
                tree.exec_on_segment(l, r, [b](auto& m) {m.chmax = b;},
                                     meta::proceed_chmax(b), meta::stop_chmax(b));
                for(int i = l; i < r; ++i) expected[i] = std::max(expected[i], b);
            }
            int64_t sum = 0;
            tree.exec_on_segment(l, r, [&](auto& m) {sum += m.sum;});
            assert(sum == std::accumulate(expected.begin() + l, expected.begin() + r, int64_t(0)));
        }
    }
}
