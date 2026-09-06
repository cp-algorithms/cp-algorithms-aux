#ifndef CP_ALGO_MATH_POLY_SQRT_HPP
#define CP_ALGO_MATH_POLY_SQRT_HPP
#include "inv.hpp"
#include "../../number_theory/discrete_sqrt.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math {
    // A square root modulo x^n, or nullopt if none exists.
    template<typename T>
    std::optional<poly_t<T>> sqrt(poly_t<T> p, size_t n) {
        if(n == 0) {return poly_t<T>{};}
        p.mod_xk_inplace(n);
        if(p.is_zero()) {return p;}
        size_t shift = p.trailing_xk();
        if(shift % 2) {return std::nullopt;}
        if(shift) {
            p.div_xk_inplace(shift);
            auto ans = sqrt(std::move(p), n - shift / 2);
            if(ans) {ans->mul_xk_inplace(shift / 2);}
            return ans;
        }
        auto c = math::sqrt(p[0]);
        if(!c) {return std::nullopt;}
        poly_t<T> ans = *c;
        for(size_t m = 1; m < n; m *= 2) {
            ans -= (ans - p.mod_xk(2 * m) * inv(ans, 2 * m)).mod_xk(2 * m) / 2;
        }
        ans.mod_xk_inplace(n);
        return ans;
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_POLY_SQRT_HPP
