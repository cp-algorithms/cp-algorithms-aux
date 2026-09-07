#include "cp-algo/math/poly/recurrence.hpp"
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
    auto monic = [](P p) {return p.is_zero() ? p : p / p.lead();};
    for(size_t n: {0, 1, 2, 15, 63, 64, 65, 127, 128, 129, 257}) {
        for(int trial = 0; trial < 5; trial++) {
            auto common = random(1 + rng() % 13), a = random(n), b = random(1 + rng() % 150);
            a *= common; b *= common;
            auto x = a, y = b;
            while(!y.is_zero()) {
                auto r = poly::impl::divmod_slow(std::move(x), y)[1];
                x = std::move(y); y = std::move(r);
            }
            auto want = monic(x);
            assert(monic(gcd(a, b)) == want);
            assert(monic(gcd(b, a)) == want);
            auto inverse = inv_mod(a, b);
            assert(bool(inverse) == (want.deg() == 0));
            if(inverse && b.deg() > 0) {assert((a * *inverse) % b == P(1));}
        }
    }
    for(int d: {1, 2, 3, 7, 31, 32, 33, 65}) {
        auto q = random(d + 1); q.a[0] = 1;
        auto seq = inv(q, 2*d + 5);
        auto r = min_rec(seq, 2*d + 5);
        assert(r.deg() <= d);
        for(int start = 0; start + r.deg() < 2*d + 5; start++) {
            T sum = 0;
            for(int j = 0; j <= r.deg(); j++) {sum += r[j] * seq[start+j];}
            assert(sum == T(0));
        }
    }
    assert(gcd(P{}, P{}).is_zero());
    assert(min_rec(P{}, 100) == P(1));
    for(size_t n = 1; n <= 65; n++) {
        for(size_t at = 0; at < n; at++) {
            assert(monic(min_rec(P::xk(at), n)) == P::xk(at+1));
        }
    }
}
int main() {
    check<modint<998244353>>();
    check<modint<1000000007>>();
    std::cout << "Polynomial GCD, modular inverse, and minimal recurrence properties passed\n";
}
