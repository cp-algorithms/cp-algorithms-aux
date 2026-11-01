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
        assert(apply_diff_even(p, h) == expected);
        assert(apply_diff_even(p, h) == apply_diff(p, h.x2()));
        assert(apply_diff_even(p, poly_t<T>{}).is_zero());
        auto integrated = integr_apply_diff_even(p, h);
        assert(integrated[0] == T(0));
        assert(deriv(integrated) == expected);
        assert(integr_apply_diff_even(p, poly_t<T>{}).is_zero());
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
