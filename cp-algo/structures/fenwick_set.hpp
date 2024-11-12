#ifndef CP_ALGO_STRUCTURES_FENWICK_SET_HPP
#define CP_ALGO_STRUCTURES_FENWICK_SET_HPP
#include "fenwick.hpp"
#include <immintrin.h>
#include <cstdint>
#include <bitset>
namespace cp_algo::structures {
    // fenwick-based set for [0, maxc)
    // Requires GCC target("popcnt,bmi2")
    template<typename Uint>
    constexpr size_t width = sizeof(Uint) * 8;
    template<size_t maxc, typename Uint>
    using popcount_array = std::array<int, maxc / width<Uint> + 1>;
    template<size_t maxc, typename Uint = uint64_t>
    struct fenwick_set: fenwick<int, popcount_array<maxc, Uint>> {
        using Base = fenwick<int, popcount_array<maxc, Uint>>;
        static constexpr size_t word = width<Uint>;
        size_t sz = 0;
        std::array<Uint, maxc / word + 1> bits;

        void flip_bit(size_t x) {
            bits[x / word] ^= 1ULL << (x % word);
        }
        bool present(size_t x) const {
            return (bits[x / word] >> (x % word)) & 1;
        }

        fenwick_set(): Base(popcount_array<maxc, Uint>{}) {}
        fenwick_set(auto &&range): fenwick_set() {
            for(auto x: range) {
                Base::data[x / word + 1] += 1;
                if(!present(x)) {
                    sz++;
                    flip_bit(x);
                }
            }
            Base::to_prefix_sums();
        }
        void insert(size_t x) {
            if(present(x)) return;
            flip_bit(x);
            sz++;
            Base::add(x / word, 1);
        }
        void erase(size_t x) {
            if(!present(x)) return;
            flip_bit(x);
            sz--;
            Base::add(x / word, -1);
        }
        static size_t order_of_bit(Uint x, size_t k) {
            return k ? std::popcount(x << (word - k)) : 0;
        }
        size_t order_of_key(size_t x) const {
            return Base::prefix_sum(x / word) + order_of_bit(bits[x / word], x % word);
        }
        static size_t kth_set_bit(Uint x, size_t k) {
            return std::countr_zero(_pdep_u64(1ULL << k, x));
        }
        size_t find_by_order(size_t order) const {
            if(order >= sz) {
                return -1;
            }
            auto [x, remainder] = Base::prefix_lower_bound(order + 1);
            return x * word + kth_set_bit(bits[x], remainder - 1);
        }
        size_t lower_bound(size_t x) const {
            if(present(x)) {return x;}
            auto order = order_of_key(x);
            return order < sz ? find_by_order(order) : -1;
        }
        size_t pre_upper_bound(size_t x) const {
            if(present(x)) {return x;}
            auto order = order_of_key(x);
            return order ? find_by_order(order - 1) : -1;
        }
    };
}
#endif // CP_ALGO_STRUCTURES_FENWICK_SET_HPP
