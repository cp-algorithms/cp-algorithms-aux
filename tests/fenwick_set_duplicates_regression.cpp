#include <cassert>
#include <bits/stdc++.h>
#ifdef TEST_MINIFIED
#include "cp-algo/min/structures/fenwick_set.hpp"
#else
#include "cp-algo/structures/fenwick_set.hpp"
#endif
int main() {
    std::mt19937 rng(7328);
    for(int test = 0; test < 1000; ++test) {
        std::vector<int> input = {1, 1, 70};
        for(int i = 0, n = rng() % 300; i < n; ++i)
            input.push_back(rng() % 128);
        cp_algo::structures::fenwick_set<128> actual(input);
        std::set<int> expected(input.begin(), input.end());
        assert(actual.sz == expected.size());
        for(int x = 0; x < 128; ++x) {
            assert(actual.order_of_key(x) == size_t(std::distance(expected.begin(), expected.lower_bound(x))));
        }
        size_t k = 0;
        for(int x : expected)
            assert(actual.find_by_order(k++) == size_t(x));
        assert(actual.find_by_order(k) == size_t(-1));
    }
}
