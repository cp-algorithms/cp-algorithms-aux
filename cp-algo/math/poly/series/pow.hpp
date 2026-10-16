#ifndef CP_ALGO_MATH_POLY_SERIES_POW_HPP
#define CP_ALGO_MATH_POLY_SERIES_POW_HPP
#include "exp.hpp"
#include "log.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math {
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
#endif // CP_ALGO_MATH_POLY_SERIES_POW_HPP
