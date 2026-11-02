#ifndef CP_ALGO_MATH_POLY_TRANSFORM_HPP
#define CP_ALGO_MATH_POLY_TRANSFORM_HPP
#include "calculus.hpp"
#include "series/inv.hpp"
#include <ranges>
#include <vector>
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
    template<std::ranges::view V, typename Stride>
    struct dilated_view {
        V coefficients;
        Stride stride;
    };
    // Represent p(x^stride) without inserting zero coefficients.
    template<typename P>
    auto dilate(P&& p, size_t stride) {
        assert(stride > 0);
        return dilated_view{std::views::all(std::forward<P>(p).a), stride};
    }
    template<size_t Stride, typename P>
    auto dilate(P&& p) {
        static_assert(Stride > 0);
        return dilated_view{std::views::all(std::forward<P>(p).a), std::integral_constant<size_t, Stride>{}};
    }
namespace poly::impl {
    // Own the correlations; defer conversion from factorial basis.
    template<typename T, typename Stride = size_t>
    struct diff_result {
        std::vector<poly_t<T>> parts;
        Stride stride;

        poly_t<T> eval(int offset = 0) const {
            size_t n = 0;
            for(size_t j = 0; j < parts.size(); ++j) {
                if(!parts[j].is_zero()) {n = std::max(n, stride * parts[j].deg() + j + 1);}
            }
            if(n == 0) {return {};}
            typename poly_t<T>::Vector result(n + offset);
            for(size_t j = 0; j < parts.size(); ++j) {
                for(int i = 0; i <= parts[j].deg(); ++i) {
                    int k = int(stride * i + j) + offset;
                    result[k] = parts[j][i] * rfact<T>(k);
                }
            }
            return result;
        }
        operator poly_t<T>() const& {return eval();}
        operator poly_t<T>() && {
            if(stride == 1 && !parts.empty()) {return borel(std::move(parts[0]));}
            return eval();
        }
    };
}
    template<typename T, typename Stride>
    poly_t<T> integr(poly::impl::diff_result<T, Stride> const& p) {return p.eval(1);}

    // g(D)p; materialize as poly_t, or integrate before materializing.
    // Input and output factorials must be invertible and within the cache.
    template<typename T, std::ranges::view V, typename Stride>
    auto apply_diff(poly_t<T> const& p, dilated_view<V, Stride> g) {
        size_t count = std::min(size_t(g.stride), p.a.size());
        poly::impl::diff_result<T, Stride> result{{}, g.stride};
        result.parts.resize(count);
        typename poly_t<T>::Vector coefficients(std::ranges::begin(g.coefficients), std::ranges::end(g.coefficients));
        poly_t<T> h(std::move(coefficients));
        for(size_t j = 0; j < count; ++j) {
            auto &part = result.parts[j];
            part.a.resize((p.a.size() - 1 - j) / g.stride + 1);
            for(size_t i = 0, k = j; i < part.a.size(); ++i, k += g.stride) {part.a[i] = p[k] * fact<T>(int(k));}
            part.normalize();
            part = forward_corr(std::move(part), j + 1 == count ? std::move(h) : h);
        }
        return result;
    }
    template<typename T>
    auto apply_diff(poly_t<T> p, poly_t<T> g) {
        poly::impl::diff_result<T> result{{}, 1};
        result.parts.push_back(forward_corr(invborel(std::move(p)), std::move(g)));
        return result;
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
