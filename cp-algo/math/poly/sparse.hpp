#ifndef CP_ALGO_MATH_POLY_SPARSE_HPP
#define CP_ALGO_MATH_POLY_SPARSE_HPP
#include "base.hpp"
#include "../combinatorics.hpp"
#include "../../number_theory/discrete_sqrt.hpp"
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
namespace cp_algo::math {
    // Inverse modulo x^n in O(n * number of nonzero terms), with p[0] != 0.
    template<typename T>
    poly_t<T> inv_sparse(poly_t<T> const& p, size_t n) {
        if(n == 0) {return {};}
        assert(p[0] != T(0));
        auto terms = poly::impl::sparse_terms(p, n);
        T a0inv = T(1) / p[0];
        for(auto &[j, a]: terms) {a *= -a0inv;}
        typename poly_t<T>::Vector q(n);
        q[0] = a0inv;
        for(size_t i = 1; i < n; i++) {
            for(auto [j, a]: terms) {
                if(j > i) {break;}
                q[i] += a * q[i - j];
            }
        }
        return q;
    }
    // Exponential modulo x^n in O(n * number of nonzero terms), with p[0] = 0.
    template<typename T>
    poly_t<T> exp_sparse(poly_t<T> const& p, size_t n) {
        if(n == 0) {return {};}
        assert(p[0] == T(0));
        auto terms = poly::impl::sparse_terms(p, n);
        for(auto &[j, a]: terms) {a *= T(j);}
        typename poly_t<T>::Vector q(n);
        q[0] = 1;
        for(size_t i = 1; i < n; i++) {
            for(auto [j, a]: terms) {
                if(j > i) {break;}
                q[i] += a * q[i - j];
            }
            q[i] *= small_inv<T>(i);
        }
        return q;
    }
    // Logarithm modulo x^n in O(n * number of nonzero terms), with p[0] = 1.
    template<typename T>
    poly_t<T> log_sparse(poly_t<T> const& p, size_t n) {
        if(n == 0) {return {};}
        assert(p[0] == T(1));
        auto terms = poly::impl::sparse_terms(p, n);
        typename poly_t<T>::Vector q(n);
        // Store x * log(p)' first, so the convolution has an unweighted kernel.
        for(size_t i = 1; i < n; i++) {
            q[i] = T(i) * p[int(i)];
            for(auto [j, a]: terms) {
                if(j > i) {break;}
                q[i] -= a * q[i - j];
            }
        }
        for(size_t i = 1; i < n; i++) {q[i] *= small_inv<T>(i);}
        return q;
    }
    // Nonnegative integer power modulo x^n in O(n * number of nonzero terms).
    template<typename T>
    poly_t<T> pow_sparse(poly_t<T> const& p, int64_t k, size_t n) {
        assert(k >= 0);
        if(n == 0) {return {};}
        if(k == 0) {return T(1);}
        if(p.is_zero()) {return {};}
        size_t shift = p.trailing_xk();
        if(shift && uint64_t(k) > (n - 1) / shift) {return {};}
        size_t offset = shift * uint64_t(k);
        auto q = poly::impl::pow_sparse_unit(p, T(k), n - offset, shift);
        q *= bpow(p.a[shift], k);
        q.mul_xk_inplace(offset);
        return q;
    }
    // Square root modulo x^n in O(n * number of nonzero terms), or nullopt.
    template<typename T>
    std::optional<poly_t<T>> sqrt_sparse(poly_t<T> const& p, size_t n) {
        if(n == 0 || p.is_zero()) {return poly_t<T>{};}
        size_t shift = p.trailing_xk();
        if(shift >= n) {return poly_t<T>{};}
        if(shift % 2) {return std::nullopt;}
        auto c = math::sqrt(p.a[shift]);
        if(!c) {return std::nullopt;}
        auto q = poly::impl::pow_sparse_unit(p, T(1) / T(2), n - shift, shift);
        q *= *c;
        q.mul_xk_inplace(shift / 2);
        return q;
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_POLY_SPARSE_HPP
