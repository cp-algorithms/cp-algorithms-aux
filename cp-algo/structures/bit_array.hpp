#ifndef CP_ALGO_STRUCTURES_BIT_ARRAY_HPP
#define CP_ALGO_STRUCTURES_BIT_ARRAY_HPP
#include "cp-algo/util/bit.hpp"
namespace cp_algo::structures {
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
#endif // CP_ALGO_STRUCTURES_BIT_ARRAY_HPP
