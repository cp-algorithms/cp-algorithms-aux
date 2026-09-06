#ifndef CP_ALGO_MATH_POLY_INV_HPP
#define CP_ALGO_MATH_POLY_INV_HPP
#include "base.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::poly::impl {
    template<typename poly>
    poly& inv_inplace(poly& p, size_t n) {
        using poly_t = std::decay_t<poly>;
        using base = poly_t::base;
        if(n == 0) {
            p.a.clear();
            return p;
        }
        assert(p[0] != base(0));
        if(n == 1) {
            return p = base(1) / p[0];
        }
        // P(x) = q0(x^2) + x q1(x^2).
        auto [q0, q1] = p.bisect(n);

        size_t N = fft::com_size((n + 1) / 2, (n + 1) / 2);

        auto q0f = fft::dft<base>(q0.a, N);
        auto q1f = fft::dft<base>(q1.a, N);

        // Q(x)*Q(-x) = Q0(x^2)^2 - x^2 Q1(x^2)^2
        auto qq = poly_t(q0f * q0f) - poly_t(q1f * q1f).mul_xk_inplace(1);

        inv_inplace(qq, (n + 1) / 2);
        auto qqf = fft::dft<base>(qq.a, N);

        typename poly::Vector A, B;
        A.resize(((n + 1) / 2 + fft::flen - 1) / fft::flen * fft::flen);
        B.resize(((n + 1) / 2 + fft::flen - 1) / fft::flen * fft::flen);
        q0f.mul(qqf, A, (n + 1) / 2);
        q1f.mul_inplace(qqf, B, (n + 1) / 2);
        p.a.resize(n + 1);
        for(size_t i = 0; i < n; i += 2) {
            p.a[i] = A[i / 2];
            p.a[i + 1] = -B[i / 2];
        }
        p.a.pop_back();
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
#endif // CP_ALGO_MATH_POLY_INV_HPP
