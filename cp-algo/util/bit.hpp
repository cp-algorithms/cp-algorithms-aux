#ifndef CP_ALGO_UTIL_BIT_HPP
#define CP_ALGO_UTIL_BIT_HPP
#include <immintrin.h>
#include <cstdint>
#include <array>
#include <bit>
namespace cp_algo {
    template<typename Uint>
    constexpr size_t bit_width = sizeof(Uint) * 8;

    size_t order_of_bit(auto x, size_t k) {
        return k ? std::popcount(x << (bit_width<decltype(x)> - k)) : 0;
    }
    // Requires GCC target("popcnt,bmi2")
    size_t kth_set_bit(uint64_t x, size_t k) {
        return std::countr_zero(_pdep_u64(1ULL << k, x));
    }
}
#endif // CP_ALGO_UTIL_BIT_HPP
