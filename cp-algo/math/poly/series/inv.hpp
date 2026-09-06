#ifndef CP_ALGO_MATH_POLY_SERIES_INV_HPP
#define CP_ALGO_MATH_POLY_SERIES_INV_HPP
#include "../base.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::poly::impl {
    template<typename poly>
    poly& inv_inplace(poly& p, size_t n) {
        using base = poly::base;
        if(n == 0) {
            p.a.clear();
            return p;
        }
        assert(p[0] != base(0));
        if(n < magic) {
            typename poly::Vector q(n);
            q[0] = base(1) / p[0];
            for(size_t i = 1; i < n; i++) {
                for(size_t j = 1; j <= std::min(i, p.a.size() - 1); j++) {
                    q[i] -= p.a[j] * q[i - j];
                }
                q[i] *= q[0];
            }
            return p = std::move(q);
        }
        size_t m = std::bit_floor(size_t(magic - 1));
        auto q = p.mod_xk(m);
        inv_inplace(q, m);
        for(; m < n; m *= 2) {
            size_t k = std::min(2 * m, n);
            typename poly::Vector error((k + fft::flen - 1) / fft::flen * fft::flen);
            auto Q = fft::dft<base>(q.a, m);
            {
                auto P = fft::dft<base>(p.a | std::views::take(k), m);
                // Wrapping modulo x^(2m) + factor^(2m) only changes the discarded low half.
                P.mul(Q, error, k);
            }
            auto E = fft::dft<base>(error | std::views::drop(m) | std::views::take(k - m), m);
            Q.mul_inplace(E, error, k - m);
            q.a.resize(k);
            for(size_t i = m; i < k; i++) {q.a[i] = -error[i - m];}
        }
        p = std::move(q);
        p.normalize();
        return p;
    }
}
namespace cp_algo::math {
    // Inverse modulo x^n; the constant coefficient must be invertible.
    template<typename T>
    poly_t<T> inv(poly_t<T> p, size_t n) {
        poly::impl::inv_inplace(p, n);
        return p;
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_POLY_SERIES_INV_HPP
