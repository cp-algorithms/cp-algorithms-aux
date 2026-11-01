#ifndef CP_ALGO_MATH_POLY_TRANSFORM_HPP
#define CP_ALGO_MATH_POLY_TRANSFORM_HPP
#include "calculus.hpp"
#include "series/inv.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math {
    // Multiply coefficient k by c^k, i.e. substitute cx for x.
    template<typename T>
    poly_t<T> mulx(poly_t<T> p, std::type_identity_t<T> c) {
        T cur = 1;
        for(auto &x: p.a) {x *= cur; cur *= c;}
        p.normalize();
        return p;
    }
    // Multiply coefficient k by c^{k(k+1)/2} (chirp factors).
    template<typename T>
    poly_t<T> mulx_sq(poly_t<T> p, std::type_identity_t<T> c) {
        T cur = 1, total = 1;
        for(auto &x: p.a) {x *= total; cur *= c; total *= cur;}
        p.normalize();
        return p;
    }
    template<typename T>
    poly_t<T> invborel(poly_t<T> p) { // a[k] *= k!
        for(int i = 0; i <= p.deg(); i++) {p.a[i] *= fact<T>(i);}
        return p;
    }
    template<typename T>
    poly_t<T> borel(poly_t<T> p) { // a[k] /= k!
        for(int i = 0; i <= p.deg(); i++) {p.a[i] *= rfact<T>(i);}
        return p;
    }
    template<typename T>
    poly_t<T> expx(size_t n) {return borel(poly_t<T>::ones(n));}

    template<typename T>
    poly_t<T> log1px(size_t n) {
        typename poly_t<T>::Vector a(n);
        for(size_t i = 1; i < n; i++) {a[i] = (i & 1 ? small_inv<T>(i) : -small_inv<T>(i));}
        return a;
    }
    template<typename T>
    poly_t<T> log1mx(size_t n) {
        return n ? -integr(poly_t<T>::ones(n - 1)) : poly_t<T>{};
    }
    // Cross-correlation, with the reversed second argument.
    template<typename T>
    poly_t<T> corr(poly_t<T> a, poly_t<T> b) {
        b.reverse();
        a *= b;
        return a;
    }
    // Coefficient k is sum_i a[i+k] b[i].
    template<typename T>
    poly_t<T> forward_corr(poly_t<T> a, poly_t<T> b) {
        if(b.is_zero()) {return {};}
        int d = b.deg();
        return corr(std::move(a), std::move(b)).div_xk(d);
    }
    // Dot product of b with every full-length substring of a.
    template<typename T>
    poly_t<T> inner_corr(poly_t<T> a, poly_t<T> b) {
        if(b.is_zero() || a.deg() < b.deg()) {return {};}
        int d = b.deg(), n = a.deg() - d + 1;
        size_t N = std::max(fft::flen, std::bit_ceil(a.a.size()));
        std::ranges::reverse(b.a);
        a.a.resize(N);
        b.a.resize(N);
        fft::cyclic_mul(a.a, b.a, N);
        a.substr_inplace(d, n);
        return a;
    }
    template<typename T>
    poly_t<T> apply_diff(poly_t<T> p, poly_t<T> g) { // g(D) p(x)
        return borel(forward_corr(invborel(std::move(p)), std::move(g)));
    }
    template<typename T>
    poly_t<T> shift(poly_t<T> p, std::type_identity_t<T> c) { // p(x+c)
        auto g = mulx(expx<T>(p.deg() + 1), c);
        return apply_diff(std::move(p), std::move(g));
    }
    // q(0) = 0 and q(x+1) - q(x) = p(x).
    template<typename T>
    poly_t<T> prefix_sum(poly_t<T> p) {
        if(p.is_zero()) {return {};}
        size_t n = p.a.size(), m = (n + 1) / 2;
        typename poly_t<T>::Vector c(m), s(m), result(n + 1);
        T half = T(1) / T(2), quarter = half * half, scale = 1;
        // t/(exp(t)-1) = (t/2)*coth(t/2) - t/2: only the linear term is odd.
        for(size_t j = 0; j < m; ++j) {
            c[j] = scale * rfact<T>(2 * j);
            s[j] = scale * rfact<T>(2 * j + 1);
            scale *= quarter;
        }
        auto kernel = poly_t<T>(std::move(c));
        kernel.mul_truncate(inv(poly_t<T>(std::move(s)), m), m);
        p = invborel(std::move(p));
        for(size_t parity = 0; parity < 2; ++parity) {
            typename poly_t<T>::Vector a((n + 1 - parity) / 2);
            for(size_t j = 0; j < a.size(); ++j) {a[j] = p.a[2 * j + parity];}
            auto corr = forward_corr(poly_t<T>(std::move(a)), kernel);
            for(size_t k = parity; k < n; k += 2) {
                result[k + 1] = (corr[k / 2] - half * p[k + 1]) * rfact<T>(k + 1);
            }
        }
        return result;
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_POLY_TRANSFORM_HPP
