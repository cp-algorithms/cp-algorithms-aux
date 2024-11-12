#ifndef CP_ALGO_STRUCTURES_FENWICK_SET_HPP
#define CP_ALGO_STRUCTURES_FENWICK_SET_HPP
#include "fenwick.hpp"
#include <bitset>
namespace cp_algo::structures {
    // fenwick-based set for [0, maxc)
    template<size_t maxc>
    struct fenwick_set: fenwick<int, std::array<int, maxc+1>> {
        using Base = fenwick<int, std::array<int, maxc+1>>;
        size_t sz = 0;
        std::bitset<maxc> present;
        fenwick_set(): Base(std::array<int, maxc+1>()) {}
        fenwick_set(auto &&range): fenwick_set() {
            for(auto x: range) {
                Base::data[x + 1] = 1;
                sz += !present[x];
                present[x] = 1;
            }
            Base::to_prefix_sums();
        }
        void insert(size_t x) {
            if(present[x]) return;
            present[x] = 1;
            sz++;
            Base::add(x, 1);
        }
        void erase(size_t x) {
            if(!present[x]) return;
            present[x] = 0;
            sz--;
            Base::add(x, -1);
        }
        size_t order_of_key(size_t x) const {
            return Base::prefix_sum(x);
        }
        size_t find_by_order(size_t order) const {
            return order < sz ? Base::prefix_lower_bound(order + 1) : -1;
        }
        size_t lower_bound(size_t x) const {
            if(present[x]) {return x;}
            auto order = order_of_key(x);
            return order < sz ? find_by_order(order) : -1;
        }
        size_t pre_upper_bound(size_t x) const {
            if(present[x]) {return x;}
            auto order = order_of_key(x);
            return order ? find_by_order(order - 1) : -1;
        }
    };
}
#endif // CP_ALGO_STRUCTURES_FENWICK_SET_HPP
