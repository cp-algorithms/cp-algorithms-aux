#include "cp-algo/math/poly/recurrence.hpp"
#include <random>
#include <iostream>
using namespace cp_algo::math;
template<typename T> size_t bm_degree(std::vector<T> const& a) {
    std::vector<T> c{1}, b{1};
    size_t len = 0, shift = 1;
    T previous = 1;
    for(size_t n = 0; n < a.size(); n++) {
        T error = a[n];
        for(size_t i = 1; i <= len; i++) {if(i < c.size()) {error += c[i]*a[n-i];}}
        if(error == T(0)) {shift++; continue;}
        auto old = c;
        T ratio = error / previous;
        c.resize(std::max(c.size(), b.size()+shift));
        for(size_t i = 0; i < b.size(); i++) {c[i+shift] -= ratio*b[i];}
        if(2*len <= n) {len = n+1-len; b = std::move(old); previous = error; shift = 1;}
        else {shift++;}
    }
    return len;
}
template<typename T> void check() {
    using P = poly_t<T>;
    std::mt19937 rng(193);
    auto random = [&](size_t n) {
        typename P::Vector a(n);
        for(auto &x: a) {x = rng() % T::mod();}
        return P(std::move(a));
    };
    auto monic = [](P p) {return p.is_zero() ? p : p / p.lead();};
    auto recurrence = [&](std::vector<T> const& a) {
        auto r = min_rec(P(typename P::Vector(a.begin(),a.end())), a.size());
        assert(r.deg() == int(bm_degree(a)));
        for(size_t i = 0; i+size_t(r.deg()) < a.size(); i++) {
            T value = 0;
            for(int j = 0; j <= r.deg(); j++) {value += r[j]*a[i+j];}
            assert(value == T(0));
        }
    };
    for(size_t n = 0; n <= 11; n++) {
        for(size_t mask = 0; mask < (size_t(1)<<n); mask++) {
            std::vector<T> a(n);
            for(size_t i = 0; i < n; i++) {a[i] = (mask>>i)&1;}
            recurrence(a);
        }
    }
    for(size_t n: {31,32,33,63,64,65,127,128,129,255,256,257,513}) {
        for(int rep = 0; rep < 8; rep++) {
            std::vector<T> a(n);
            for(auto &x:a) {x = rng()%T::mod();}
            recurrence(a);
        }
    }

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
    auto mul = [](P const& a, P const& b) {
        typename P::Vector c(a.a.size()+b.a.size());
        for(size_t i=0;i<a.a.size();i++)for(size_t j=0;j<b.a.size();j++){c[i+j]+=a.a[i]*b.a[j];}
        return P(std::move(c));
    };
    for(size_t n: {1,2,3,31,63,64,65,127,129}) {
        for(size_t m: {0,1,2,31,63,64,65,129}) {
            for(int monic = 0; monic < 2; monic++) {
                auto q=random(n);q.a.back()=monic?T(1):T(17);
                auto quotient=random(m), rem=random(n-1);
                auto dividend=mul(quotient,q)+rem;
                auto [d,r]=divmod(dividend,q);
                assert(d==quotient && r==rem);
            }
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
