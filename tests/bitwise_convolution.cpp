#include <bits/stdc++.h>
#include "cp-algo/number_theory/modint.hpp"
#include "cp-algo/math/xor_convolution.hpp"
#include "cp-algo/math/and_convolution.hpp"
using namespace cp_algo::math;
template<class T> void check() {
    std::mt19937 rng(81378);
    for(int rep = 0; rep < 60; rep++) {
        size_t n = 1u << (rep % 9);
        std::vector<T> a(n), b(n), want_xor(n), want_and(n);
        for(auto &x: a) x = rng() % T::mod();
        for(auto &x: b) x = rng() % T::mod();
        for(size_t i = 0; i < n; i++) for(size_t j = 0; j < n; j++) {
            want_xor[i ^ j] += a[i] * b[j];
            want_and[i & j] += a[i] * b[j];
        }
        assert(xor_convolution(a, b) == want_xor);
        assert(and_convolution(a, b) == want_and);
    }
}
int main() {
    check<modint<998244353>>();
    check<modint<1000000007>>();
    dynamic_modint<>::with_mod(998244353, [] {check<dynamic_modint<>>();});
    std::cout << "Bitwise convolutions passed naive checks at both recursion parities\n";
}
