#ifndef CP_ALGO_MATH_POLY_CALCULUS_HPP
#define CP_ALGO_MATH_POLY_CALCULUS_HPP
#include "base.hpp"
#include "../factorials.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math {
    // k-th derivative; consumes temporaries and moved polynomials.
    template<typename T>
    poly_t<T> deriv(poly_t<T> p, int k = 1) {
        assert(k >= 0);
        if(k > p.deg()) {return k == 0 ? p : poly_t<T>{};}
        for(int i = k; i <= p.deg(); i++) {
            p.a[i - k] = fact<T>(i) * rfact<T>(i - k) * p.a[i];
        }
        p.a.resize(p.a.size() - k);
        p.normalize();
        return p;
    }
    // Antiderivative with zero constant coefficient.
    template<typename T>
    poly_t<T> integr(poly_t<T> p) {
        if(p.is_zero()) {return p;}
        p.a.push_back(0);
        for(int i = p.deg() - 1; i >= 0; i--) {
            p.a[i + 1] = p.a[i] * small_inv<T>(i + 1);
        }
        p.a[0] = 0;
        return p;
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_POLY_CALCULUS_HPP
