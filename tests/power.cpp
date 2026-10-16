#include "cp-algo/math/poly.hpp"
#include "cp-algo/linalg/matrix.hpp"
#include <random>
#include <iostream>
using namespace cp_algo::math;
struct nondefault {
    uint64_t value;
    nondefault() = delete;
    explicit nondefault(uint64_t v): value(v) {}
};
template<int W> void counts(uint64_t n) {
    size_t binary = 0, window = 0;
    auto a = bpow(nondefault(1), n, nondefault(0), [&](auto a, auto b) {
        binary++; return nondefault(a.value + b.value);
    });
    auto b = bpow<W>(nondefault(1), n, nondefault(0), [&](auto a, auto b) {
        window++; return nondefault(a.value + b.value);
    });
    assert(a.value == n && b.value == n && window <= binary);
}
int main() {
    std::mt19937_64 rng(42);
    for(uint64_t n = 0; n < 65536; n++) {counts<2>(n); counts<3>(n); counts<4>(n); counts<5>(n); counts<6>(n);}
    for(int i = 0; i < 10000; i++) {
        auto n = rng(); counts<3>(n); counts<4>(n);
        using T = modint<998244353>;
        T x = rng();
        assert(bpow<3>(x, n, T(1)) == bpow(x, n));
    }
    counts<1>(UINT64_MAX); counts<2>(UINT64_MAX); counts<3>(UINT64_MAX);
    counts<4>(UINT64_MAX); counts<5>(UINT64_MAX); counts<6>(UINT64_MAX);
    for(int i = 0; i < 64; i++) {counts<3>(uint64_t(1) << i);}
    using T = modint<998244353>; using P = poly_t<T>;
    for(int n: {1, 2, 3, 16, 65, 129}) {
        P::Vector a(n), b(n + 1);
        for(auto &x: a) {x = rng() % T::mod();}
        for(auto &x: b) {x = rng() % T::mod();}
        b.back() = 1;
        for(P md: {P(b), P::xk(n) - P(T(1)), P::xk(n)}) {
            for(int64_t k: {0LL, 1LL, 2LL, 3LL, 15LL, 31LL, 63LL, 1023LL, 1000000000000000000LL}) {
                auto want = bpow(P(a) % md, k, P(T(1)), [&](auto const& a, auto const& b) {return a * b % md;});
                assert(powmod(P(a), k, md) == want);
            }
        }
        assert(powmod(P{}, 5, P(b)).is_zero());
    }
    using U = modint<int64_t(998244353)>; using M = cp_algo::linalg::matrix<U>;
    for(int n: {0, 1, 2, 5}) {
        M a(n, n), want = M::eye(n);
        for(int i = 0; i < n; i++) {for(int j = 0; j < n; j++) {a[i][j] = rng() % U::mod();}}
        for(uint64_t k: {uint64_t(1000000000000000000), uint64_t(1) << 63, UINT64_MAX}) {
            auto actual = a.pow(k), expected = bpow(a, k, M::eye(n));
            for(int i = 0; i < n; i++) {for(int j = 0; j < n; j++) {assert(actual[i][j] == expected[i][j]);}}
        }
        for(int k = 0; k <= 40; k++) {
            auto actual = a.pow(k);
            for(int i = 0; i < n; i++) {for(int j = 0; j < n; j++) {assert(actual[i][j] == want[i][j]);}}
            want = want * a;
        }
    }
    std::cout << "Window counts, 64-bit exponents, nondefault monoids, matrices and polynomial powers passed\n";
}
