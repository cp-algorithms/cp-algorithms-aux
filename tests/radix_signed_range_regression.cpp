#include <cassert>
#include <bits/stdc++.h>
#ifdef TEST_MINIFIED
#include "cp-algo/min/util/sort.hpp"
#include "cp-algo/min/util/compress_coords.hpp"
#else
#include "cp-algo/util/sort.hpp"
#include "cp-algo/util/compress_coords.hpp"
#endif
template <class T> void check(std::vector<T> values) {
    auto expected = values;
    std::sort(expected.begin(), expected.end());
    cp_algo::radix_sort(values);
    assert(values == expected);
}
int main() {
    check<int>({2000000000, -2000000000});
    check<int>({INT_MAX, 0, INT_MIN, -1, INT_MAX});
    check<unsigned>({UINT_MAX, 0, UINT_MAX / 2, 1});
    __int128 high = __int128(1) << 100;
    check<__int128>({high + 2, high, high + 1});
    std::mt19937 rng(8731);
    for(int test = 0; test < 1000; ++test) {
        std::vector<int> values(1 + rng() % 256);
        for(auto &x : values)
            x = int(rng());
        check(values);
    }
    std::vector<int> values = {INT_MAX, INT_MIN, 0, INT_MAX};
    std::vector<std::reference_wrapper<int>> refs(values.begin(), values.end());
    auto original = cp_algo::compress_coords(refs);
    assert((values == std::vector<int>{2, 0, 1, 2}));
    assert((std::vector<int>(original.begin(), original.end()) == std::vector<int>{INT_MIN, 0, INT_MAX}));
}
