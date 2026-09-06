#ifndef CP_ALGO_MATH_POLY_CHIRPZ_HPP
#define CP_ALGO_MATH_POLY_CHIRPZ_HPP
#include "transform.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math {
    // p(1), p(z), ..., p(z^{n-1}); needs convolution length n + 2 deg(p).
    template<typename T>
    poly_t<T> chirpz(poly_t<T> p, std::type_identity_t<T> z, int n) {
        assert(n >= 0);
        if(n == 0 || p.is_zero()) {return {};}
        if(z == T(0)) {
            typename poly_t<T>::Vector ans(n, p[0]);
            ans[0] = p.eval(1);
            return ans;
        }
        auto b = mulx_sq(poly_t<T>::ones(n + p.deg()), z);
        auto a = mulx_sq(std::move(p), z.inv());
        return mulx_sq(inner_corr(std::move(b), std::move(a)), z.inv());
    }
namespace poly::impl {
    // res[i] = product_{1 <= j <= i} 1/(1-z^j).
    template<typename T>
    big_vector<T> geometric_products_inv(T z, int n) {
        big_vector<T> res(n, 1), zk(n, 1);
        for(int i = 1; i < n; i++) {
            zk[i] = zk[i - 1] * z;
            res[i] = res[i - 1] * (T(1) - zk[i]);
        }
        res.back() = res.back().inv();
        for(int i = n - 2; i >= 0; i--) {res[i] = (T(1) - zk[i + 1]) * res[i + 1];}
        return res;
    }
    // product_{0 <= j < n} (1-z^j x).
    template<typename T>
    poly_t<T> geometric_product(T z, int n) {
        if(n == 0) {return T(1);}
        if(n == 1) {return typename poly_t<T>::Vector{1, -1};}
        auto p = geometric_product(z, n / 2);
        p *= mulx(p, bpow(z, n / 2));
        if(n % 2) {p *= poly_t<T>(typename poly_t<T>::Vector{1, -bpow(z, n - 1)});}
        return p;
    }
}
    // Interpolate from values at 1,z,...,z^{n-1}; nonzero nodes must be distinct.
    template<typename T>
    poly_t<T> chirpz_inverse(poly_t<T> p, std::type_identity_t<T> z, int n) {
        assert(n >= 0);
        if(n == 0 || p.is_zero()) {return {};}
        if(z == T(0)) {
            return n == 1 ? p.mod_xk(1) : poly_t<T>(typename poly_t<T>::Vector{p[1], p[0] - p[1]});
        }
        p.a.resize(n);
        auto pos = poly::impl::geometric_products_inv(z, n);
        auto neg = poly::impl::geometric_products_inv(z.inv(), n);
        T zn = bpow(z, n - 1).inv(), znk = 1;
        for(int i = 0; i < n; i++) {
            p.a[i] *= znk * neg[i] * pos[n - 1 - i];
            znk *= zn;
        }
        p = chirpz(std::move(p), z, n);
        p *= poly::impl::geometric_product(z, n);
        p.mod_xk_inplace(n).reverse(n);
        return p;
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_POLY_CHIRPZ_HPP
