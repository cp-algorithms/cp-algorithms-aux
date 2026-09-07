#include "cp-algo/math/poly/sparse.hpp"
#include "cp-algo/math/poly/series.hpp"
#include "cp-algo/math/poly/series/sqrt.hpp"
#include <random>
#include <iostream>
using namespace cp_algo::math;
template<typename T>
void check() {
    using P = poly_t<T>;
    std::mt19937 rng(42);
    for(int n: {0, 1, 2, 3, 31, 32, 63, 64, 65, 127, 129, 257}) {
        for(int trial = 0; trial < 15; trial++) {
            typename P::Vector a(n + 17);
            for(int j = 0; j < 6; j++) {a[rng() % a.size()] = rng() % T::mod();}
            P p(a);
            for(int64_t k: {int64_t(0), int64_t(1), int64_t(2), int64_t(7), int64_t(1000000000000000000)}) {
                assert(pow_sparse(p, k, n) == pow(p, k, n));
            }
            auto q = sqrt_sparse(p, n), want = sqrt(p, n);
            assert(bool(q) == bool(want));
            if(q) {assert((*q * *q).mod_xk(n) == p.mod_xk(n));}
            a[0] = 17;
            p = a;
            assert(inv_sparse(p, n) == inv(p, n));
            a[0] = 1;
            p = a;
            assert(log_sparse(p, n) == log(p, n));
            a[0] = 0;
            p = a;
            assert(exp_sparse(p, n) == exp(p, n));
        }
    }
    for(int n: {0, 1, 2, 7, 64, 129}) {
        assert(pow_sparse(P{}, 0, n) == (n ? P(1) : P{}));
        assert(pow_sparse(P{}, 13, n).is_zero());
        assert(sqrt_sparse(P{}, n)->is_zero());
        assert(inv_sparse(P(3), n) == (n ? P(T(1) / T(3)) : P{}));
        assert(log_sparse(P(1), n).is_zero());
        assert(exp_sparse(P{}, n) == (n ? P(1) : P{}));
        for(int shift = 0; shift < n + 3; shift++) {
            P p = P({4, 0, 9}).mul_xk(shift);
            auto q = sqrt_sparse(p, n);
            assert(bool(q) == (shift >= n || shift % 2 == 0));
            if(q) {assert((*q * *q).mod_xk(n) == p.mod_xk(n));}
            assert(pow_sparse(p, 1000000000000000000, n) == pow(p, 1000000000000000000, n));
        }
    }
}
template<typename T> void accumulation_boundaries() {
    using P = poly_t<T>;
    constexpr size_t n = 80;
    auto mul = [](P const& a, P const& b) {
        typename P::Vector c(n);
        for(size_t i = 0; i < a.a.size(); i++) {
            for(size_t j = 0; j < b.a.size() && i+j < n; j++) {c[i+j] += a.a[i] * b.a[j];}
        }
        return P(std::move(c));
    };
    for(size_t terms: {0, 1, 7, 8, 9, 15, 16, 17, 31}) {
        typename P::Vector a(terms+1, T::mod()-1);
        a[0] = 1;
        P p(a);
        assert(mul(p, inv_sparse(p, n)) == P(1));
        auto l = log_sparse(p, n);
        assert(mul(p, deriv(l)).mod_xk(n-1) == deriv(p).mod_xk(n-1));
        for(int64_t k: {int64_t(2), int64_t(7), int64_t(1000000000000000000)}) {
            auto q = pow_sparse(p, k, n);
            assert(q[0] == T(1));
            assert(mul(p, deriv(q)).mod_xk(n-1) == (mul(deriv(p), q) * T(k)).mod_xk(n-1));
        }
        auto root = sqrt_sparse(p, n);
        assert(root && mul(*root, *root) == p.mod_xk(n));
        p.a[0] = 0;
        auto e = exp_sparse(p, n);
        assert(e[0] == T(1));
        assert(mul(deriv(p), e).mod_xk(n-1) == deriv(e));
    }
}
int main() {
    check<modint<998244353>>();
    check<modint<1000000007>>();
    accumulation_boundaries<modint<998244353>>();
    accumulation_boundaries<modint<2013265921>>();
    accumulation_boundaries<modint<998244353LL>>();
    dynamic_modint<>::with_mod(998244353, [] {accumulation_boundaries<dynamic_modint<>>();});
    std::cout << "Sparse polynomial properties passed\n";
}
