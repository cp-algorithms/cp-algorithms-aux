#ifndef CP_ALGO_STRUCTURES_FENWICK_SET_HPP
#define CP_ALGO_STRUCTURES_FENWICK_SET_HPP
#include "fenwick.hpp"
#include "cp-algo/bit.hpp"
namespace cp_algo::structures {
    // fenwick-based set for [0, maxc)
    template<size_t maxc, typename Uint = uint64_t>
    struct fenwick_set: fenwick<int, popcount_array<maxc, Uint>> {
        using Base = fenwick<int, popcount_array<maxc, Uint>>;
        static constexpr size_t word = bit_width<Uint>;
        size_t sz = 0;
        bit_array<maxc, Uint> bits;

        fenwick_set(): Base(popcount_array<maxc, Uint>{}) {}
        fenwick_set(auto &&range): fenwick_set() {
            for(auto x: range) {
                Base::data[x / word + 1] += 1;
                if(!bits.test(x)) {
                    sz++;
                    bits.flip(x);
                }
            }
            Base::to_prefix_sums();
        }
        void insert(size_t x) {
            if(bits.test(x)) return;
            bits.flip(x);
            sz++;
            Base::add(x / word, 1);
        }
        void erase(size_t x) {
            if(!bits.test(x)) return;
            bits.flip(x);
            sz--;
            Base::add(x / word, -1);
        }
        size_t order_of_key(size_t x) const {
            return Base::prefix_sum(x / word) + order_of_bit(bits.word(x / word), x % word);
        }
        size_t find_by_order(size_t order) const {
            if(order >= sz) {
                return -1;
            }
            auto [x, remainder] = Base::prefix_lower_bound(order + 1);
            return x * word + kth_set_bit(bits.word(x), remainder - 1);
        }
        size_t lower_bound(size_t x) const {
            if(bits.test(x)) {return x;}
            auto order = order_of_key(x);
            return order < sz ? find_by_order(order) : -1;
        }
        size_t pre_upper_bound(size_t x) const {
            if(bits.test(x)) {return x;}
            auto order = order_of_key(x);
            return order ? find_by_order(order - 1) : -1;
        }
    };
}
#endif // CP_ALGO_STRUCTURES_FENWICK_SET_HPP
