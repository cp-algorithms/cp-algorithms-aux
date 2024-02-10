#ifndef CP_ALGO_DATA_STRUCTURES_TREAP_METAS_REVERSE_HPP
#define CP_ALGO_DATA_STRUCTURES_TREAP_METAS_REVERSE_HPP
#include "base.hpp"
#include <algorithm>
namespace cp_algo::data_structures::treap::metas {
    struct reverse_meta: base_meta {
        int val;
        bool reverse = false;
        int64_t sum = val;

        reverse_meta(int val): val(val) {}

        void pull(auto const L, auto const R) {
            sum = val + _safe_meta(L, sum) + _safe_meta(R, sum);
        }
        void push(auto &L, auto &R) {
            if(reverse) {
                reverse = false;
                std::swap(L, R);
                _safe_meta(L, reverse ^= 1);
                _safe_meta(R, reverse ^= 1);
            }
        }
    };
}
#endif // CP_ALGO_DATA_STRUCTURES_TREAP_METAS_REVERSE_HPP