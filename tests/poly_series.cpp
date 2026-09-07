#include "cp-algo/math/poly.hpp"
#include <random>
#include <iostream>
using namespace cp_algo::math;
template<typename T> void check() {
    using P = poly_t<T>;
    std::mt19937 rng(193);
    auto random = [&](size_t n) {
        typename P::Vector a(n);
        for(auto &x: a) {x = rng() % T::mod();}
        return P(std::move(a));
    };
    auto multiply = [](P const& a, P const& b, size_t n) {
        typename P::Vector c(n);
        for(size_t i = 0; i < a.a.size(); i++) {
            for(size_t j = 0; j < b.a.size() && i + j < n; j++) {c[i + j] += a.a[i] * b.a[j];}
        }
        return P(std::move(c));
    };
    for(size_t n: {1, 2, 31, 32, 33, 63, 64, 65, 127, 128, 129, 255, 256, 257, 511, 512, 513, 1025, 1537}) {
        assert(log(P(T(1)), n).is_zero());
        auto f = random(n); f.a[0] = 0;
        auto e = exp(f, n);
        assert(e[0] == T(1));
        assert(deriv(e) == multiply(deriv(f), e, n - 1));
        assert(log(e, n) == f.mod_xk(n));
        auto a = random(n); a.a[0] = 17;
        assert(deriv(a, 0) == a);
        assert(deriv(a, 2) == deriv(deriv(a)));
        auto square = multiply(a, a, n);
        auto root = sqrt(square, n);
        assert(root && multiply(*root, *root, n) == square);
        for(size_t shift: {size_t(0), size_t(2), 2 * (n / 3), n}) {
            auto p = square.mul_xk(shift).mod_xk(n);
            auto r = sqrt(p, n);
            assert(r && multiply(*r, *r, n) == p);
        }
    }
    for(size_t n: {1, 2, 31, 32, 33, 63, 64, 65, 129, 257}) {
        typename P::Vector x(n), y(n);
        auto a = random(n);
        for(size_t i = 0; i < n; i++) {x[i] = T(i) * T(7) + T(123); y[i] = a.eval(x[i]);}
        assert(eval(a, x) == y);
        assert(inter(x, y) == a);
        std::fill(begin(y), end(y), T(0));
        assert(inter(x, y).is_zero());
    }
    // Include z of exact order n: 1-z^n is zero, but all interpolation nodes differ.
    for(size_t n: {1, 2, 4, 16, 32, 64, 128, 256}) {
        if((T::mod() - 1) % n) {continue;}
        T z = n == 2 ? T(-1) : bpow(T(3), (T::mod() - 1) / n), x = 1;
        auto a = random(n);
        typename P::Vector y(n);
        for(size_t i = 0; i < n; i++, x *= z) {y[i] = a.eval(x);}
        assert(chirpz_inverse(P(y), z, n) == a);
    }
    for(size_t m: {1, 2, 64, 65, 127, 128, 129, 256}) {
        for(size_t n: {size_t(0), m - 1, m, m + 1, m + 65, m + 257, 1537ul}) {
            auto a = random(n), b = random(m);
            auto [q, r] = divmod(a, b);
            assert(r.deg() < b.deg());
            assert(multiply(q, b, n) + r == a);
        }
    }
}
int main() {
    check<modint<998244353>>();
    check<modint<1000000007>>();
    std::cout << "series, interpolation and division properties passed under two primes\n";
}
