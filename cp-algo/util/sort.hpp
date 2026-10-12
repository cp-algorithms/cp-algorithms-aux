#ifndef CP_ALGO_UTIL_SORT_HPP
#define CP_ALGO_UTIL_SORT_HPP
#include "bit.hpp"
#include <algorithm>
#include <numeric>
#include <ranges>
#include <type_traits>
#include <vector>
namespace cp_algo {
    template<size_t maxc>
    void count_sort(auto &a, auto &&proj) {
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
    template<size_t maxc>
    void count_sort(auto &a) {
        count_sort<maxc>(a, std::identity{});
    }

    void radix_sort(auto &a, auto &&proj) {
        if(empty(a)) {
            return;
        }
        auto [mn, mx] = std::ranges::minmax(a, {}, proj);
        using value_type = decltype(std::invoke(proj, mx) - std::invoke(proj, mn));
        using key_type = typename std::conditional_t<
            (sizeof(value_type) > sizeof(uint64_t)),
            std::type_identity<unsigned __int128>,
            std::make_unsigned<value_type>>::type;
        auto key = [&](auto const& x) -> key_type {
            return key_type(std::invoke(proj, x)) - key_type(std::invoke(proj, mn));
        };
        with_bit_floor<1>(size(a), [&]<size_t floor>() {
            constexpr int base = std::min<size_t>(floor, 1 << 16);
            for(int64_t i = 1; i <= key(mx); i *= base) {
                count_sort<base>(a, [&](auto const& x) {
                    return key(x) / i % base;
                });
            }
        });
    }
    void radix_sort(auto &a) {
        radix_sort(a, std::identity{});
    }
}
#endif // CP_ALGO_UTIL_SORT_HPP
