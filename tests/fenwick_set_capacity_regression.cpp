#include <cassert>
#include <bits/stdc++.h>
#ifdef TEST_MINIFIED
#include "cp-algo/min/structures/fenwick_set.hpp"
#else
#include "cp-algo/structures/fenwick_set.hpp"
#endif
template <size_t N> void check() {
    for(size_t first = 0; first < N; ++first) {
        std::vector<int> values;
        for(size_t x = first; x < N; ++x)
            values.push_back(x);
        cp_algo::structures::fenwick_set<N> actual(values);
        assert(actual.sz == values.size());
        for(size_t x = 0; x < N; ++x) {
            assert(actual.order_of_key(x) == (x < first ? 0 : x - first));
        }
        for(size_t k = 0; k < values.size(); ++k)
            assert(actual.find_by_order(k) == size_t(values[k]));
        actual.erase(N - 1);
        assert(actual.sz + 1 == values.size());
        actual.insert(N - 1);
        assert(actual.find_by_order(actual.sz - 1) == N - 1);
    }
}
int main() {
    check<1>();
    check<2>();
    check<63>();
    check<64>();
    check<65>();
    check<127>();
    check<128>();
    check<129>();
    check<257>();
}
