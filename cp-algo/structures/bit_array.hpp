#ifndef CP_ALGO_STRUCTURES_BIT_ARRAY_HPP
#define CP_ALGO_STRUCTURES_BIT_ARRAY_HPP
#include "../util/bit.hpp"
#include "../util/bump_alloc.hpp"
namespace cp_algo::structures {
    template<class Cont>
    struct _bit_array {
        static constexpr size_t width = bit_width<uint64_t>;
        size_t words, n;
        alignas(32) Cont data;

        _bit_array(): words(0), n(0) {}
        _bit_array(size_t N): words((N + width - 1) / width), n(N), data() {}

        uint64_t& word(size_t x) {
            return data[x];
        }
        uint64_t word(size_t x) const {
            return data[x];
        }
        void set(size_t x) {
            word(x / width) |= 1ULL << (x % width);
        }
        void reset(size_t x) {
            word(x / width) &= ~(1ULL << (x % width));
        }
        void reset() {
            for(auto& w: data) {
                w = 0;
            }
        }
        void flip(size_t x) {
            word(x / width) ^= 1ULL << (x % width);
        }
        bool test(size_t x) const {
            return (word(x / width) >> (x % width)) & 1;
        }
        bool operator[](size_t x) const {
            return test(x);
        }
        size_t size() const {
            return n;
        }
    };

    template<int N>
    struct bit_array: _bit_array<std::array<uint64_t, (N + 63) / 64>> {
        using Base = _bit_array<std::array<uint64_t, (N + 63) / 64>>;
        using Base::Base, Base::words, Base::data;
        bit_array(): Base(N) {
            data.fill(0);
        }
    };
    struct dynamic_bit_array: _bit_array<std::vector<uint64_t>> {
        using Base = _bit_array<std::vector<uint64_t>>;
        using Base::Base, Base::words;
        dynamic_bit_array(size_t N): Base(N) {
            data.resize(words);
        }
    };

}
#endif // CP_ALGO_STRUCTURES_BIT_ARRAY_HPP
