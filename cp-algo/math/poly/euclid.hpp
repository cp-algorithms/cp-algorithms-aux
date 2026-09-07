#ifndef CP_ALGO_MATH_POLY_EUCLID_HPP
#define CP_ALGO_MATH_POLY_EUCLID_HPP
#include "impl/euclid.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math {
    // GCD, without normalizing the leading coefficient.
    template<typename T>
    poly_t<T> gcd(poly_t<T> a, poly_t<T> b) {
        poly::impl::full_gcd<false, false>(a, b);
        return a;
    }
    // Inverse modulo q, or nullopt when p and q are not coprime.
    template<typename T>
    std::optional<poly_t<T>> inv_mod(poly_t<T> p, poly_t<T> q) {
        return poly::impl::inv_mod(std::move(p), std::move(q));
    }
    template<typename T>
    T resultant(poly_t<T> a, poly_t<T> b) {
        T res = 1;
        while(!b.is_zero()) {
            if(b.deg() == 0) {return res * bpow(b.lead(), a.deg());}
            int d = a.deg();
            a %= b;
            res *= bpow(b.lead(), d - a.deg()) * T((b.deg() & a.deg() & 1) ? -1 : 1);
            std::swap(a, b);
        }
        return T(0);
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_POLY_EUCLID_HPP
