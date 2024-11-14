#ifndef CP_ALGO_UTIL_SORT_HPP
#define CP_ALGO_UTIL_SORT_HPP
#include <algorithm>
#include <numeric>
#include <ranges>
#include <vector>
namespace cp_algo {
    void count_sort(auto &a, size_t maxc, auto &&proj = std::identity{}) {
        std::vector<int> cnt(maxc);
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
        int base = std::bit_floor(size(a));
        auto mx = std::ranges::max(a);
        for(int64_t i = 1; i <= mx; i *= base) {
            count_sort(a, base, [&](auto x) {
                return x / i % base;
            });
        }
    }
}
#endif // CP_ALGO_UTIL_SORT_HPP
