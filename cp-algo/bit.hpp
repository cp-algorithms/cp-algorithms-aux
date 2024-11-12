#ifndef CP_ALGO_BIT_HPP
#define CP_ALGO_BIT_HPP
#include <immintrin.h>
#include <cstdint>
#include <array>
#include <bit>
namespace cp_algo {
    template<typename Uint>
    constexpr size_t bit_width = sizeof(Uint) * 8;
    template<size_t maxc, typename Uint = uint64_t>
    using popcount_array = std::array<int, maxc / bit_width<Uint> + 1>;

    size_t order_of_bit(auto x, size_t k) {
        return k ? std::popcount(x << (bit_width<decltype(x)> - k)) : 0;
    }
    // Requires GCC target("popcnt,bmi2")
    size_t kth_set_bit(uint64_t x, size_t k) {
        return std::countr_zero(_pdep_u64(1ULL << k, x));
    }

    template<size_t N, typename Uint = uint64_t>
    struct bit_array {
        static constexpr size_t width = bit_width<Uint>;
        static constexpr size_t blocks = N / width + 1;
        std::array<Uint, blocks> data = {};

        uint64_t word(size_t x) const {
            return data[x];
        }
        void set(size_t x) {
            data[x / width] |= 1ULL << (x % width);
        }
        void flip(size_t x) {
            data[x / width] ^= 1ULL << (x % width);
        }
        bool test(size_t x) const {
            return (data[x / width] >> (x % width)) & 1;
        }
        bool operator[](size_t x) const {
            return test(x);
        }
    };
}
#endif // CP_ALGO_BIT_HPP
