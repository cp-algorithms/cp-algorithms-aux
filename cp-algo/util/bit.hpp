#ifndef CP_ALGO_UTIL_BIT_HPP
#define CP_ALGO_UTIL_BIT_HPP
#include "../util/simd.hpp"
#include <cstdint>
#include <array>
#include <bit>
namespace cp_algo {
    template<typename Uint>
    constexpr size_t bit_width = sizeof(Uint) * 8;

    // n < 64
    uint64_t mask(size_t n) {
        return (1ULL << n) - 1;
    }
    size_t order_of_bit(auto x, size_t k) {
        return std::popcount(x << ( -k % bit_width<decltype(x)>));
    }
    [[gnu::target("bmi2")]] inline size_t kth_set_bit(uint64_t x, size_t k) {
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
    void with_bit_ceil(size_t n, auto &&callback) {
        with_bit_floor(n, [&]<size_t N>() {
            if(N == n) {
                callback.template operator()<N>();
            } else {
                callback.template operator()<N << 1>();
            }
        });
    }

    [[gnu::target("avx2")]] inline uint32_t read_bits(char const* p) {
        return _mm256_movemask_epi8(__m256i(vector_cast<u8x32 const>(p[0]) + (127 - '0')));
    }
    [[gnu::target("avx2")]] inline uint64_t read_bits64(char const* p) {
        return read_bits(p) | (uint64_t(read_bits(p + 32)) << 32);
    }

    [[gnu::target("avx2")]] inline void write_bits(char *p, uint32_t bits) {
        static constexpr u8x32 shuffler = {
            0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3
        };
        auto shuffled = u8x32(_mm256_shuffle_epi8(__m256i() + bits, __m256i(shuffler)));
        static constexpr u8x32 mask = {
            1, 2, 4, 8, 16, 32, 64, 128,
            1, 2, 4, 8, 16, 32, 64, 128,
            1, 2, 4, 8, 16, 32, 64, 128,
            1, 2, 4, 8, 16, 32, 64, 128
        };
        for(int z = 0; z < 32; z++) {
            p[z] = shuffled[z] & mask[z] ? '1' : '0';
        }
    }
    [[gnu::target("avx2")]] inline void write_bits64(char *p, uint64_t bits) {
        write_bits(p, uint32_t(bits));
        write_bits(p + 32, uint32_t(bits >> 32));
    }
}
#endif // CP_ALGO_UTIL_BIT_HPP
