#include <bits/stdc++.h>
#ifdef TEST_MINIFIED
#include "cp-algo/min/util/sort.hpp"
#else
#include "cp-algo/util/sort.hpp"
#endif
template<class T>
void check(std::vector<T> values) {
    auto expected = values;
    std::sort(expected.begin(), expected.end());
    cp_algo::radix_sort(values);
    assert(values == expected);
}
int main() {
    std::vector<int64_t> values(256);
    values[0] = 1000000000000000000LL;
    check(values);
    check<int64_t>({INT64_MAX, INT64_MIN, 0, -1, 1});
    check<uint64_t>({UINT64_MAX, 0, UINT64_MAX / 2, 1});
    __int128 lowest = -(__int128(1) << 126) * 2;
    check<__int128>({-(lowest + 1), lowest, 0, -1, 1});
    check<unsigned __int128>({~static_cast<unsigned __int128>(0), 0, 1});
    std::mt19937_64 rng(8351);
    for(size_t n: {2, 3, 4, 7, 16, 63, 64, 127, 128, 255, 256, 257}) {
        for(int test = 0; test < 20; ++test) {
            std::vector<int64_t> signed_values(n);
            std::vector<uint64_t> unsigned_values(n);
            for(auto& x: signed_values) x = int64_t(rng());
            for(auto& x: unsigned_values) x = rng();
            signed_values[0] = INT64_MIN;
            signed_values[1] = INT64_MAX;
            unsigned_values[0] = 0;
            unsigned_values[1] = UINT64_MAX;
            check(signed_values);
            check(unsigned_values);
        }
    }
}
