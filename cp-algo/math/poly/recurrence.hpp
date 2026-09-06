#ifndef CP_ALGO_MATH_POLY_RECURRENCE_HPP
#define CP_ALGO_MATH_POLY_RECURRENCE_HPP
#include "euclid.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::poly::impl {
    template<typename poly>
    poly& inv_inplace(poly& q, int64_t k, size_t n) {
        if(n == 0) {q.a.clear(); return q;}
        assert(k >= 0 || uint64_t(-(k + 1)) < n);
        using poly_t = std::decay_t<poly>;
        using base = poly_t::base;
        if(k <= std::max<int64_t>(n, size(q.a))) {
            inv_inplace(q, size_t(k + int64_t(n)));
            return q.div_xk_inplace(k);
        }
        if(k % 2) {
            return inv_inplace(q, k - 1, n + 1).div_xk_inplace(1);
        }
        auto [q0, q1] = q.bisect();
        auto qq = q0 * q0 - (q1 * q1).mul_xk_inplace(1);
        inv_inplace(qq, k / 2 - q.deg() / 2, (n + 1) / 2 + q.deg() / 2);
        size_t N = fft::com_size(size(q0.a), size(qq.a));
        auto q0f = fft::dft<base>(q0.a, N);
        auto q1f = fft::dft<base>(q1.a, N);
        auto qqf = fft::dft<base>(qq.a, N);
        size_t M = q0.deg() + (n + 1) / 2;
        typename poly::Vector A, B;
        A.resize((M + fft::flen - 1) / fft::flen * fft::flen);
        B.resize((M + fft::flen - 1) / fft::flen * fft::flen);
        q0f.mul(qqf, A, M);
        q1f.mul_inplace(qqf, B, M);
        q.a.resize(n + 1);
        for(size_t i = 0; i < n; i += 2) {
            q.a[i] = A[q0.deg() + i / 2];
            q.a[i + 1] = -B[q0.deg() + i / 2];
        }
        q.a.pop_back();
        q.normalize();
        return q;
    }
}
namespace cp_algo::math {
    // Non-monic characteristic polynomial of a minimum recurrence.
    template<typename T>
    poly_t<T> min_rec(poly_t<T> const& p, size_t d) {return poly::impl::min_rec(p, d);}
    // Coefficients [x^k] through [x^{k+n-1}] of 1/p; k+n must be nonnegative.
    template<typename T>
    poly_t<T> inv(poly_t<T> p, int64_t k, size_t n) {
        poly::impl::inv_inplace(p, k, n);
        return p;
    }
    // Find [x^k] P / Q
    template<typename T>
    T kth_rec(poly_t<T> P, poly_t<T> Q, int64_t k) {
        assert(k >= 0 && Q[0] != T(0));
        while(k > Q.deg()) {
            size_t n = Q.a.size();
            auto [Q0, Q1] = Q.bisect();
            auto [P0, P1] = P.bisect();

            size_t N = fft::com_size((n + 1) / 2, (n + 1) / 2);

            auto Q0f = fft::dft<T>(Q0.a, N);
            auto Q1f = fft::dft<T>(Q1.a, N);
            auto P0f = fft::dft<T>(P0.a, N);
            auto P1f = fft::dft<T>(P1.a, N);

            Q = poly_t<T>(Q0f * Q0f) - poly_t<T>(Q1f * Q1f).mul_xk_inplace(1);
            if(k % 2) {
                P = poly_t<T>(Q0f *= P1f) - poly_t<T>(Q1f *= P0f);
            } else {
                P = poly_t<T>(Q0f *= P0f) - poly_t<T>(Q1f *= P1f).mul_xk_inplace(1);
            }
            k /= 2;
        }
        size_t n = size_t(k) + 1;
        P.mul_truncate(inv(std::move(Q), n), n);
        return P[(int)k];
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_POLY_RECURRENCE_HPP
