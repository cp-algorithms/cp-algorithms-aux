#ifndef CP_ALGO_STRUCTURES_BIT_ARRAY_HPP
#define CP_ALGO_STRUCTURES_BIT_ARRAY_HPP
#include "../util/bit.hpp"
#include "../util/bump_alloc.hpp"
#include <cassert>
namespace cp_algo::structures {
    template<typename C>
    concept Resizable = requires(C& c, std::size_t n) { c.resize(n); };

    template<class Cont>
    struct _bit_array {
        static constexpr size_t width = bit_width<uint64_t>;
        size_t words, n;
        alignas(32) Cont data;

        void resize(size_t N) {
            n = N;
            words = (n + width - 1) / width;
            if constexpr (Resizable<Cont>) {
                data.resize(words);
            } else {
                assert(std::size(data) >= words);
            }
        }

        _bit_array(): n(0), words(0), data() {}
        _bit_array(size_t N): data() {
            resize(N);
        }

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
        bit_array(): Base(N) {}
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
