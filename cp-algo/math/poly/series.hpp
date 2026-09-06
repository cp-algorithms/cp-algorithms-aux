#ifndef CP_ALGO_MATH_POLY_SERIES_HPP
#define CP_ALGO_MATH_POLY_SERIES_HPP
#include "inv.hpp"
#include "calculus.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math {
    // log(p) modulo x^n, for p[0] = 1.
    template<typename T>
    poly_t<T> log(poly_t<T> p, size_t n) {
        if(n == 0) {return {};}
        assert(p[0] == T(1));
        p.mod_xk_inplace(n);
        auto dp = deriv(p);
        poly::impl::inv_inplace(p, n);
        p.mul_truncate(dp, n - 1);
        return integr(std::move(p));
    }
    // exp(p) modulo x^n, for p[0] = 0.
    template<typename T>
    poly_t<T> exp(poly_t<T> p, size_t n) {
        if(n == 0) {return {};}
        assert(p[0] == T(0));
        p.mod_xk_inplace(n);
        if(p.is_zero()) {return T(1);}
        p.a[0] = 1;
        for(size_t m = 1; m < n; m *= 2) {
            auto c = log(p, 2 * m).div_xk(m) - p.substr(m, 2 * m);
            c.mul_truncate(p, m).mul_xk_inplace(m);
            p -= c;
        }
        p.mod_xk_inplace(n);
        return p;
    }
    namespace poly::impl {
        // O(deg(p) * n), using p q' = k p' q.
        template<typename T>
        poly_t<T> pow_slow(poly_t<T> const& p, int64_t k, size_t n) {
            typename poly_t<T>::Vector q(n);
            q[0] = bpow(p[0], k);
            auto a0inv = p[0].inv();
            for(int i = 1; i < (int)n; i++) {
                for(int j = 1; j <= std::min(p.deg(), i); j++) {
                    q[i] += p[j] * q[i - j] * (T(k) * T(j) - T(i - j));
                }
                q[i] *= small_inv<T>(i) * a0inv;
            }
            return q;
        }
    }
    // Nonnegative integer power modulo x^n.
    template<typename T>
    poly_t<T> pow(poly_t<T> p, int64_t k, size_t n) {
        assert(k >= 0);
        if(n == 0) {return {};}
        if(k == 0) {return T(1);}
        p.mod_xk_inplace(n);
        if(p.is_zero()) {return p;}
        size_t shift = p.trailing_xk();
        if(shift) {
            if(uint64_t(k) > (n - 1) / shift) {return {};}
            p.div_xk_inplace(shift);
            return pow(std::move(p), k, n - shift * k).mul_xk(shift * k);
        }
        if(std::min(p.deg(), (int)n) <= magic) {
            return poly::impl::pow_slow(p, k, n);
        }
        if(k <= magic) {
            auto t = pow(p, k / 2, n);
            t.mul_truncate(t, n);
            if(k % 2) {t.mul_truncate(p, n);}
            return t;
        }
        T c = p[0];
        p /= c;
        return bpow(c, k) * exp(log(std::move(p), n) * T(k), n);
    }

}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_POLY_SERIES_HPP
