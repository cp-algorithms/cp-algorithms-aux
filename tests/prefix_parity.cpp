#include "cp-algo/math/poly/transform.hpp"
#include <cassert>
#include <random>
using namespace cp_algo::math;
template<class T> void check() {
    std::mt19937 rng(42);
    for(int n = 0; n <= 48; ++n) {
        typename poly_t<T>::Vector a(n);
        for(auto &x : a) x = rng() % T::mod();
        poly_t<T> p(a);
        typename poly_t<T>::Vector hv(n / 2 + 2);
        for(auto &x: hv) x = rng() % T::mod();
        poly_t<T> h(hv), expected;
        auto d = p;
        for(auto c: hv) {
            expected += d * c;
            d = deriv(deriv(d));
        }
        for(size_t stride: {1, 2, 3, 4, 7, 64}) {
            auto expected = poly_t<T>{};
            auto d = p;
            for(auto c: hv) {expected += d * c; d = deriv(d, int(stride));}
            poly_t<T> actual = apply_diff(p, dilate(h, stride));
            assert(actual == expected);
            auto integrated = integr(apply_diff(p, dilate(h, stride)));
            assert(integrated[0] == T(0));
            assert(deriv(integrated) == expected);
            assert(poly_t<T>(apply_diff(p, dilate(poly_t<T>{}, stride))).is_zero());
            auto owning = dilate(poly_t<T>(hv), stride);
            assert(poly_t<T>(apply_diff(p, std::move(owning))) == expected);
        }
        auto check_static = [&]<size_t Stride>() {
            poly_t<T> a = apply_diff(p, dilate<Stride>(h));
            poly_t<T> b = apply_diff(p, dilate(h, Stride));
            assert(a == b);
            assert(integr(apply_diff(p, dilate<Stride>(poly_t<T>(hv)))) == integr(apply_diff(p, dilate(h, Stride))));
        };
        check_static.template operator()<1>();
        check_static.template operator()<2>();
        check_static.template operator()<3>();
        check_static.template operator()<64>();
        auto d1 = p;
        auto ordinary = poly_t<T>{};
        for(auto c: hv) {ordinary += d1 * c; d1 = deriv(d1);}
        assert(poly_t<T>(apply_diff(p, h)) == ordinary);
        assert(deriv(integr(apply_diff(p, h))) == ordinary);
        auto tail = h.a | std::views::drop(1);
        poly_t<T> from_view = apply_diff(p, dilated_view{tail, size_t(2)});
        assert(from_view == poly_t<T>(apply_diff(p, dilate<2>(h.div_xk(1)))));
        auto deferred = apply_diff(poly_t<T>(a), dilate(poly_t<T>(hv), 2));
        assert(poly_t<T>(deferred) == expected);
        assert(deriv(integr(deferred)) == expected);
        auto q = prefix_sum(p);
        assert(q[0] == T(0));
        for(int x = -5; x <= 55; ++x) assert(q.eval(T(x + 1)) - q.eval(T(x)) == p.eval(T(x)));
    }
    for(int n0 = 0; n0 <= 24; ++n0) for(int n1 = 0; n1 <= 24; ++n1) {
        typename poly_t<T>::Vector a(n0), b(n1);
        for(auto &x: a) x = rng() % T::mod();
        for(auto &x: b) x = rng() % T::mod();
        poly_t<T> p0(a), p1(b);
        auto p = poly_t<T>::interleave(p0, p1);
        assert(p.bisect()[0] == p0 && p.bisect()[1] == p1);
        for(int k = 0; k < 50; ++k) assert(p[k] == (k % 2 ? p1[k / 2] : p0[k / 2]));
        assert(poly_t<T>::interleave(p.bisect()[0], p.bisect()[1]) == p);
        for(int k = 0; k < 50; ++k) assert(p0.x2()[k] == (k % 2 ? T(0) : p0[k / 2]));
    }
    assert(prefix_sum(poly_t<T>{}).is_zero());
    for(size_t n = 0; n <= 24; ++n) {
        auto c = Cx<T>(n), s = Sx<T>(n), e = expx<T>(2 * n);
        auto combined = poly_t<T>::interleave(c, s);
        assert(combined == e);
        T factorial = 1;
        for(size_t j = 0; j < 2 * n; ++j) {
            if(j) factorial *= T(j);
            assert(e[j] * factorial == T(1));
        }
    }
}
int main() {check<modint<998244353>>(); check<modint<101>>();}
