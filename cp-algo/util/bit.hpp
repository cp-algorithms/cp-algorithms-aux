#ifndef CP_ALGO_UTIL_BIT_HPP
#define CP_ALGO_UTIL_BIT_HPP
#include "../util/simd.hpp"
#include <cstdint>
#include <array>
#include <bit>
namespace cp_algo {
    template<typename Uint>
    constexpr size_t bit_width = sizeof(Uint) * 8;

    size_t order_of_bit(auto x, size_t k) {
        return k ? std::popcount(x << (bit_width<decltype(x)> - k)) : 0;
    }
    [[gnu::target("bmi2")]]
    size_t kth_set_bit(uint64_t x, size_t k) {
        return std::countr_zero(_pdep_u64(1ULL << k, x));
    }
    template<int fl = 0>
    void with_bit_floor(size_t n, auto &&callback) {
        if constexpr (fl >= 63) {
            return;
        } else if (n >> (fl + 1)) {
            with_bit_floor<fl + 1>(n, callback);
        } else {
            callback.template operator()<1ULL << fl>();
        }
    }

    [[gnu::target("avx2"), gnu::always_inline]] inline uint32_t read_bits(char const* p) {
        return _mm256_movemask_epi8(__m256i(vector_cast<u8x32 const>(p[0]) + (127 - '0')));
    }
    [[gnu::always_inline]] inline uint64_t read_bits64(char const* p) {
        return read_bits(p) | (uint64_t(read_bits(p + 32)) << 32);
    }
}
#endif // CP_ALGO_UTIL_BIT_HPP
