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
    template<typename T>
    T sparse_dot(auto const& terms, typename poly_t<T>::Vector const& q, size_t i) {
        if constexpr(T::bits <= 32) {
            // Sixteen products below 2^60 have a sum strictly below 2^64.
            if(T::mod() < (1 << 30) && terms.size() <= 16) {
                uint64_t sum = 0;
                for(auto [j, a]: terms) {
                    if(j > i) {break;}
                    sum += uint64_t(a.getr()) * q[i-j].getr();
                }
                return T(sum % T::mod());
            }
        }
        T sum = 0;
        for(auto [j, a]: terms) {
            if(j > i) {break;}
            sum += a * q[i-j];
        }
        return sum;
    }
    // For q = (p / p[shift])^k, use p q' = k p' q after removing x^shift.
    template<typename T>
    poly_t<T> pow_sparse_unit(poly_t<T> const& p, T k, size_t n, size_t shift) {
        auto terms = sparse_terms(p, n, shift);
        T a0inv = T(1) / p.a[shift];
        for(auto &[j, a]: terms) {a *= a0inv;}
        std::vector<T> weights;
        for(auto [j, a]: terms) {weights.push_back((k + T(1)) * T(j));}
        typename poly_t<T>::Vector q(n);
        q[0] = 1;
        if constexpr(T::bits <= 32) {
            if(T::mod() < (1 << 30) && terms.size() >= 6 && terms.size() <= 16) {
                for(size_t t = 0; t < terms.size(); t++) {weights[t] *= terms[t].second;}
                for(size_t i = 1; i < n; i++) {
                    uint64_t plain = 0, weighted = 0;
                    for(size_t t = 0; t < terms.size(); t++) {
                        auto [j, a] = terms[t];
                        if(j > i) {break;}
                        auto v = q[i-j].getr();
                        plain += uint64_t(a.getr()) * v;
                        weighted += uint64_t(weights[t].getr()) * v;
                    }
                    // Each sum has at most sixteen products below 2^60.
                    q[i] = (T(weighted % T::mod()) - T(i) * T(plain % T::mod())) * small_inv<T>(i);
                }
                return q;
            }
        }
        for(size_t i = 1; i < n; i++) {
            T index = T(i);
            for(size_t t = 0; t < terms.size(); t++) {
                auto [j, a] = terms[t];
                if(j > i) {break;}
                q[i] += a * q[i - j] * (weights[t] - index);
            }
            q[i] *= small_inv<T>(i);
        }
        return q;
    }
}
#pragma GCC pop_options
#endif
