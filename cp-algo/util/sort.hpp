#ifndef CP_ALGO_UTIL_SORT_HPP
#define CP_ALGO_UTIL_SORT_HPP
#include "bit.hpp"
#include <algorithm>
#include <numeric>
#include <ranges>
#include <vector>
namespace cp_algo {
    template<size_t maxc>
    void count_sort(auto &a, auto &&proj = std::identity{}) {
        std::array<int, maxc> cnt = {};
        for(auto &x: a) {
            cnt[proj(x)]++;
        }
        std::partial_sum(begin(cnt), end(cnt), begin(cnt));
        auto res = a;
        for(auto const& it: a | std::views::reverse) {
            res[--cnt[proj(it)]] = it;
        }
        a = std::move(res);
    }

    void radix_sort(auto &a) {
        if(empty(a)) {
            return;
        }
        auto mx = std::ranges::max(a);
        with_bit_floor(size(a), [&]<size_t floor>() {
            constexpr int base = std::min<size_t>(floor, 1 << 16);
            for(int64_t i = 1; i <= mx; i *= base) {
                count_sort<base>(a, [&](auto x) {
                    return x / i % base;
                });
            }
        });
    }
}
#endif // CP_ALGO_UTIL_SORT_HPP
