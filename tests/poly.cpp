#include "cp-algo/math/poly.hpp"
#include "cp-algo/math/fps.hpp"
#include "cp-algo/math/laurent.hpp"
#include <random>
#include <iostream>
using namespace cp_algo::math;
using T = modint<998244353>;
using P = poly_t<T>;
std::mt19937 rng(42);
P random_poly(int n) {
    P::Vector a(n);
    for(auto &x: a) {x = rng() % 1000;}
    return a;
}
P naive_mul(P const& a, P const& b, size_t n) {
    P::Vector c(n);
    for(int i = 0; i <= a.deg(); i++) {
        for(int j = 0; j <= b.deg() && size_t(i + j) < n; j++) {c[i + j] += a[i] * b[j];}
    }
    return c;
}
int main() {
    // Check Newton doubling at odd lengths and either side of a power of two.
    for(int n: {31, 32, 33, 63, 64, 65, 127, 128, 129, 255, 256, 257, 513, 1025}) {
        auto p = random_poly(n);
        for(auto &x: p.a) {x = rng() % T::mod();}
        p.a[0] = 17;
        assert(naive_mul(p, inv(p, n), n) == P(1));
    }
    for(int n: {0, 1, 2, 3, 7, 63, 64, 65, 127, 128, 129, 257, 513}) {
        P a = random_poly(n), b = random_poly(n);
        auto c = naive_mul(a, b, 2 * n);
        assert(a * b == c);
        auto square = a;
        square *= square;
        assert(square == naive_mul(a, a, 2 * n));
        square = a;
        square.mul_truncate(square, n);
        assert(square == naive_mul(a, a, n));
        fps<T> fa(a), fb(b);
        auto fc = fa * fb;
        for(int k = 0; k <= 2 * n; k++) {assert(fc[k] == c[k]);}
        assert(fc.prefix(2 * n) == c);
        assert((fa + fb).prefix(n) == a + b);
        assert(integr(P()).is_zero());
        if(!n) {continue;}
        a.a[0] = 1;
        auto ai = inv(a, n);
        assert(naive_mul(a, ai, n) == P(1));
        assert(inv(fps<T>(a)).prefix(n) == ai);
        assert(log(fps<T>(a)).prefix(n) == log(a, n));
        auto exponent = a; exponent.a[0] = 0;
        auto e = exp(exponent, n);
        assert(log(e, n) == exponent.mod_xk(n));
        assert(exp(fps<T>(exponent)).prefix(n) == e);
        auto root = sqrt(naive_mul(a, a, n), n);
        assert(root && naive_mul(*root, *root, n) == naive_mul(a, a, n));
        auto [q, r] = divmod(c, a);
        assert(q * a + r == c && r.deg() < a.deg());
        auto same = a; same /= same; assert(same == P(1));
        same = a; same %= same; assert(same.is_zero());
        for(int k: {0, 1, 2, 5, 65}) {
            P want(1);
            for(int j = 0; j < k; j++) {want = naive_mul(want, a, n);}
            assert(pow(a, k, n) == want);
        }
        assert(inv(a, 0).is_zero());
        assert(log(a, 0).is_zero());
        assert(exp(exponent, 0).is_zero());
        assert(pow(a, 0, 0).is_zero());
        assert(sqrt(a, 0)->is_zero());
    }
    assert(eval(P(1), P::Vector{}).empty());
    assert(inter(P::Vector{}, P::Vector{}).is_zero());
    assert(kth_rec(P(1), P({1, -1, -1}), 20) == T(10946));
    for(int k = -3; k < 20; k++) {
        auto q = inv(P({1, -1}), k, 25);
        for(int j = 0; j < 25; j++) {assert(q[j] == T(k + j >= 0));}
    }
    for(int n = 1; n < 20; n++) {
        auto a = random_poly(n), b = random_poly(4);
        b.a[0] = 0;
        P want, bk(1);
        for(int j = 0; j <= a.deg(); j++) {want += bk * a[j]; bk = naive_mul(bk, b, n);}
        assert(compose(a, b, n) == want);
        assert(compose_large(a, b, n) == want);
        P::Vector x(n), y(n);
        for(int j = 0; j < n; j++) {x[j] = j; y[j] = a.eval(j);}
        assert(inter(x, y) == a);
        auto coef = to_newton(a, x);
        want = P(); bk = P(1);
        for(int j = 0; j < n; j++) {want += bk * coef[j]; bk *= P({-T(j), 1});}
        assert(want == a);
    }
    for(int n: {2, 5, 10}) {
        auto a = random_poly(25), b = random_poly(4);
        P want, bk(1);
        for(int j = 0; j <= a.deg(); j++) {want += bk * a[j]; bk = naive_mul(bk, b, n);}
        assert(compose(a, b, n) == want);
        assert(compose_large(a, b, n) == want);
    }
    for(int n: {63, 64, 65, 129, 257}) {
        auto a = random_poly(n), b = random_poly(4);
        b.a[0] = 0;
        P want, bk(1);
        for(int j = 0; j <= a.deg(); j++) {want += bk * a[j]; bk = naive_mul(bk, b, n);}
        assert(compose(a, b, n) == want);
        assert(compose_large(a, b, n) == want);
    }
    int generated = 0;
    fps<T> fib([&](size_t n, auto const& prev) {++generated; return n < 2 ? T(1) : prev[n - 1] + prev[n - 2];});
    auto copy = fib;
    assert(fib[20] == T(10946));
    assert(copy[10] == T(89) && generated == 21);
    laurent<T> l{fps<T>(P({2, 3, 4})), -2};
    assert(l[-3] == T(0) && l[-2] == T(2) && l[0] == T(4));
    auto one = l * inv(l);
    assert(one[0] == T(1));
    for(int i = 1; i < 100; i++) {assert(one[i] == T(0));}
    auto s = l + shift(l, 3);
    for(int i = -4; i < 8; i++) {assert(s[i] == l[i] + l[i - 3]);}
    assert(deriv(l)[-3] == T(-4));
    bool rejected = false;
    try {auto unused = integr(l)[0]; (void)unused;} catch(std::domain_error const&) {rejected = true;}
    assert(rejected);
    rejected = false;
    try {auto unused = inv(fps<T>())[0]; (void)unused;} catch(std::domain_error const&) {rejected = true;}
    assert(rejected);
    fps<T> cycle;
    cycle = fps<T>([&](size_t n, auto const&) {return cycle[n + 1];});
    rejected = false;
    try {auto unused = cycle[0]; (void)unused;} catch(std::logic_error const&) {rejected = true;}
    assert(rejected);
    P movable({1, 2, 3});
    auto allocation = movable.a.data();
    auto moved = std::move(movable) * T(2);
    assert(moved.a.data() == allocation && moved == P({2, 4, 6}));
    std::cout << "polynomial, FPS and Laurent properties passed\n";
}
