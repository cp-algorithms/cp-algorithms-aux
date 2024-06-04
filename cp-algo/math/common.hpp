#ifndef CP_ALGO_MATH_COMMON_HPP
#define CP_ALGO_MATH_COMMON_HPP
#include <functional>
#include <cstdint>
namespace cp_algo::math {
#ifdef CP_ALGO_MAXN
    const int maxn = CP_ALGO_MAXN;
#else
    const int maxn = 1 << 19;
#endif
    const int magic = 128; // threshold for sizes to run the naive algo

    auto bpow(auto const& x, int64_t n, auto const& one, auto op) {
        if(n == 0) {
            return one;
        } else {
            auto t = bpow(x, n / 2, one, op);
            t = op(t, t);
            if(n % 2) {
                t = op(t, x);
            }
            return t;
        }
    }
    auto bpow(auto x, int64_t n, auto ans) {
        return bpow(x, n, ans, std::multiplies{});
    }
    template<typename T>
    T bpow(T const& x, int64_t n) {
        return bpow(x, n, T(1));
    }
}
#endif // CP_ALGO_MATH_COMMON_HPP
