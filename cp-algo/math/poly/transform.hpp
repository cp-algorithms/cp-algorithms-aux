#ifndef CP_ALGO_MATH_POLY_TRANSFORM_HPP
#define CP_ALGO_MATH_POLY_TRANSFORM_HPP
#include "calculus.hpp"
#include "coefficient_view.hpp"
#include "series/inv.hpp"
#include <ranges>
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
namespace poly::impl {
    template<typename T, size_t Step, size_t Offset>
    poly_t<T> factorial_series(size_t n) {
        typename poly_t<T>::Vector a(n);
        for(size_t i = 0; i < n; ++i) {a[i] = rfact<T>(int(Step * i + Offset));}
        return a;
    }
}
    template<typename T>
    poly_t<T> expx(size_t n) {return poly::impl::factorial_series<T, 1, 0>(n);}

    // C(x) = sum x^j/(2j)!, so cosh(x) = C(x^2).
    template<typename T>
    poly_t<T> Cx(size_t n) {return poly::impl::factorial_series<T, 2, 0>(n);}

    // S(x) = sum x^j/(2j+1)!, so sinh(x) = x*S(x^2).
    template<typename T>
    poly_t<T> Sx(size_t n) {return poly::impl::factorial_series<T, 2, 1>(n);}

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
namespace poly::impl {
    // Defer conversion from factorial basis so integration can share it.
    template<typename T>
    struct diff_result {
        poly_t<T> coefficients;
        operator poly_t<T>() const& {return borel(coefficients);}
        operator poly_t<T>() && {return borel(std::move(coefficients));}
    };
}
    // In factorial basis, integration is multiplication by x.
    template<typename T>
    poly_t<T> integr(poly::impl::diff_result<T> p) {
        return borel(std::move(p.coefficients).mul_xk(1));
    }

    // g(D)p; materialize as poly_t, or integrate before materializing.
    // Input and output factorials must be invertible and within the cache.
    template<typename T, std::ranges::input_range R>
    auto apply_diff(poly_t<T> const& p, R&& g) {
        auto [values, offset, stride] = coefficient_layout(g);
        poly_t<T> h(std::ranges::to<typename poly_t<T>::Vector>(values));
        if(offset >= p.a.size()) {return poly::impl::diff_result<T>{};}
        if(offset == 0 && stride == 1) {return poly::impl::diff_result<T>{forward_corr(invborel(p), std::move(h))};}
        size_t n = p.a.size() - offset;
        typename poly_t<T>::Vector result;
        result.reserve(n + 1); // Leave room for integration.
        result.resize(n);
        size_t count = std::min(size_t(stride), n);
        for(size_t j = 0; j < count; ++j) {
            poly_t<T> part;
            part.a.resize((n - 1 - j) / stride + 1);
            for(size_t i = 0, k = j; i < part.a.size(); ++i, k += stride) {part.a[i] = p[k + offset] * fact<T>(int(k + offset));}
            part.normalize();
            part = forward_corr(std::move(part), j + 1 == count ? std::move(h) : h);
            for(size_t i = 0, k = j; i < part.a.size(); ++i, k += stride) {result[k] = part[i];}
        }
        return poly::impl::diff_result<T>{std::move(result)};
    }
    template<typename T>
    auto apply_diff(poly_t<T> p, poly_t<T> g) {
        return poly::impl::diff_result<T>{forward_corr(invborel(std::move(p)), std::move(g))};
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
        size_t m = (p.a.size() + 1) / 2;
        T half = T(1) / T(2);
        // exp(D) shifts by one, so (exp(D)-1)q = p, D = d/dx.
        // A(t) = (exp(t)-1)/t has A(0)=1, hence q' = A(D)^(-1)p.
        // A(t)^(-1) = H(t^2)-t/2, where H(z) = C(z/4)/S(z/4).
        auto H = Cx<T>(m);
        H.mul_truncate(inv(Sx<T>(m), m), m);
        H = mulx(std::move(H), half * half);
        auto q = integr(apply_diff(p, dilate<2>(std::move(H))));
        // Subtract (p(x)-p(0))/2 to keep q(0)=0.
        for(size_t k = 1; k < p.a.size(); ++k) {q.a[k] -= half * p[k];}
        return q;
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_POLY_TRANSFORM_HPP
