#define _GLIBCXX_ASSERTIONS
#include <bits/stdc++.h>
#include "cp-algo/number_theory/modint.hpp"
#include "cp-algo/math/subset_convolution.hpp"
using namespace cp_algo;
using namespace cp_algo::math;

template<class T>
big_vector<T> naive(big_vector<T> const& a, big_vector<T> const& b) {
    big_vector<T> c(a.size());
    for(size_t s = 0; s < a.size(); s++) {
        for(size_t t = s;; t = (t - 1) & s) {
            c[s] += a[t] * b[s ^ t];
            if(!t) break;
        }
    }
    return c;
}

template<class T> void check() {
    std::mt19937 rng(87234);
    for(int rep = 0; rep < 80; rep++) {
        size_t n = size_t(1) << (rep % 8);
        big_vector<T> a(n), b(n);
        for(auto &x: a) x = rng() % T::mod();
        for(auto &x: b) x = rng() % T::mod();
        b[0] = rep + 1;
        auto c = naive(a, b);
        assert(subset_div<T>(c, b) == a);
        // The multiplication kernel supports the usual primes below 2^30.
        if(T::mod() < (1LL << 30)) {
            assert(subset_convolution<T>(a, b) == c);
            a[0] = 0;
            auto e = subset_exp<T>(a);
            assert(subset_log<T>(e) == a);
            if(n <= 32) {
                big_vector<T> f{3, 2, 7, 4}, powers(n), composed(n), projection(f.size());
                powers[0] = 1;
                for(size_t k = 0; k < f.size(); k++) {
                    for(size_t i = 0; i < n; i++) {
                        composed[i] += f[k] * powers[i];
                        projection[k] += b[i] * powers[i];
                    }
                    powers = naive(powers, a);
                }
                assert(subset_compose<T>(f, a) == composed);
                assert(subset_power_projection<T>(a, b, f.size()) == projection);
            }
        }
    }
}
int main() {
    check<modint<998244353>>();
    check<modint<1000000007>>();
    check<modint<2147483647>>();
    dynamic_modint<>::with_mod(998244353, [] {check<dynamic_modint<>>();});
    dynamic_modint<>::with_mod(1000000007, [] {check<dynamic_modint<>>();});
    using T = modint<998244353>;
    big_vector<T> coefficients{1, 2, 3, 4}, constant{2}, weight{7};
    assert(subset_compose<T>(coefficients, constant) == big_vector<T>{49});
    auto projected = subset_power_projection<T>(constant, weight, 6);
    assert(projected == (big_vector<T>{7, 14, 28, 56, 112, 224}));
    big_vector<T> f(1 << 20), g(f.size());
    for(size_t i = 0; i < f.size(); i++) {
        auto n = std::popcount(i);
        f[i] = bpow(T(-2), n);
        g[i] = bpow(T(-1), n);
    }
    assert(subset_div<T>(f, g) == g);
    std::cout << "Subset operations passed naive, changing-divisor/modulus and full-rank checks\n";
}
