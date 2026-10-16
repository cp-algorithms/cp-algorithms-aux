#ifndef CP_ALGO_MATH_COMMON_HPP
#define CP_ALGO_MATH_COMMON_HPP
#include <functional>
#include <cstdint>
#include <cassert>
#include <bit>
#include <vector>
#include <algorithm>
namespace cp_algo::math {
#ifdef CP_ALGO_MAXN
    const int maxn = CP_ALGO_MAXN;
#else
    const int maxn = 1 << 19;
#endif
    const int magic = 64; // threshold for sizes to run the naive algo

    // Nonnegative 64-bit exponents, with an associative operation and its identity.
    // Windows >1 precompute odd powers only when that saves operations.
    template<int window = 1>
    auto bpow(auto const& x, auto n, auto const& one, auto op) {
        static_assert(window >= 1 && window <= 6);
        if constexpr(window > 1) {
            if(n == 0) {return one;}
            int bits = std::bit_width(uint64_t(n));
            auto low_bit = [&](int high) {
                int low = std::max(0, high - window + 1);
                while(!((n >> low) & 1)) {low++;}
                return low;
            };
            int first = low_bit(bits - 1);
            int cost = (1 << (window - 1)) + first;
            for(int j = first - 1; j >= 0;) {
                if(!((n >> j) & 1)) {j--;}
                else {cost++; j = low_bit(j) - 1;}
            }
            // Do not pay for the table when binary powering uses fewer operations.
            if(cost >= bits + std::popcount(uint64_t(n)) - 2) {return bpow<1>(x, n, one, op);}
            using T = std::decay_t<decltype(x)>;
            std::vector<T> odd;
            odd.reserve(1 << (window - 1));
            odd.push_back(x);
            auto square = op(x, x);
            while(odd.size() < size_t(1 << (window - 1))) {odd.push_back(op(odd.back(), square));}
            auto ans = odd[(n >> first) / 2];
            for(int j = first - 1; j >= 0;) {
                if(!((n >> j) & 1)) {ans = op(ans, ans); j--;}
                else {
                    int low = low_bit(j), length = j - low + 1;
                    auto digit = (n >> low) & ((1u << length) - 1);
                    for(int i = 0; i < length; i++) {ans = op(ans, ans);}
                    ans = op(ans, odd[digit / 2]);
                    j = low - 1;
                }
            }
            return ans;
        } else {
            if (n == 0) {
                return one;
            }
            auto ans = x;
            for(int j = std::bit_width<uint64_t>(n) - 2; ~j; j--) {
                ans = op(ans, ans);
                if((n >> j) & 1) {
                    ans = op(ans, x);
                }
            }
            return ans;
        }
    }
    template<int window = 1>
    auto bpow(auto x, auto n, auto ans) {
        return bpow<window>(x, n, ans, std::multiplies{});
    }
    template<typename T>
    T bpow(T const& x, auto n) {
        return bpow(x, n, T(1));
    }
    inline constexpr auto inv2(auto x) {
        assert(x % 2);
        std::make_unsigned_t<decltype(x)> y = 1;
        while(y * x != 1) {
            y *= 2 - x * y;
        }
        return y;
    }
}
#endif // CP_ALGO_MATH_COMMON_HPP
