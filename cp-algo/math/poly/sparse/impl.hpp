#ifndef CP_ALGO_MATH_POLY_SPARSE_IMPL_HPP
#define CP_ALGO_MATH_POLY_SPARSE_IMPL_HPP
#include "../base.hpp"
#include "../../combinatorics.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::poly::impl {
    template<typename T>
    auto sparse_terms(poly_t<T> const& p, size_t n, size_t shift = 0) {
        std::vector<std::pair<size_t, T>> terms;
        for(size_t i = shift + 1; i < p.a.size() && i - shift < n; i++) {
            if(p.a[i] != T(0)) {terms.emplace_back(i - shift, p.a[i]);}
        }
        return terms;
    }
    // For q = (p / p[shift])^k, use p q' = k p' q after removing x^shift.
    template<typename T>
    poly_t<T> pow_sparse_unit(poly_t<T> const& p, T k, size_t n, size_t shift) {
        auto terms = sparse_terms(p, n, shift);
        T a0inv = T(1) / p.a[shift];
        for(auto &[j, a]: terms) {a *= a0inv;}
        typename poly_t<T>::Vector q(n);
        q[0] = 1;
        for(size_t i = 1; i < n; i++) {
            for(auto [j, a]: terms) {
                if(j > i) {break;}
                q[i] += a * q[i - j] * ((k + T(1)) * T(j) - T(i));
            }
            q[i] *= small_inv<T>(i);
        }
        return q;
    }
}
#pragma GCC pop_options
#endif
