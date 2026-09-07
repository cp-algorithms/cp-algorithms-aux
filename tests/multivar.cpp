#include "cp-algo/math/multivar.hpp"
#include <random>
#include <iostream>
using namespace cp_algo;
using namespace cp_algo::math;

template<typename T>
big_vector<T> naive(big_vector<size_t> const& dims, big_vector<T> const& a, big_vector<T> const& b) {
    big_vector<T> result(a.size());
    for(size_t i = 0; i < a.size(); i++) {
        for(size_t j = 0; i + j < a.size(); j++) {
            size_t x = i, y = j;
            bool carry = false;
            for(auto n: dims) {
                carry |= x % n + y % n >= n;
                x /= n; y /= n;
            }
            if(!carry) {result[i+j] += a[i] * b[j];}
        }
    }
    return result;
}

template<typename T> void check() {
    std::mt19937 rng(59321);
    std::vector<big_vector<size_t>> shapes{{}, {1}, {2}, {3}, {4}, {1, 2, 1, 3},
        {2, 2, 2, 2, 2, 2}, {2, 2, 2, 2, 2, 2, 2}, {3, 3, 3, 3},
        {3, 2, 3, 2, 3}, {2, 3, 2, 3, 2}, {5, 7}, {4, 3, 2}, {31, 3}, {65, 2}};
    for(int rep = 0; rep < 80; rep++) {
        big_vector<size_t> dims(rng() % 7);
        for(auto &n: dims) {n = 1 + rng() % 3;}
        shapes.push_back(dims);
    }
    for(auto const& dims: shapes) {
        fft::multivar<T> a(dims), b(dims);
        for(auto &x: a.data) {x = rng() % T::mod();}
        for(auto &x: b.data) {x = rng() % T::mod();}
        auto original = a.data, rhs = b.data;
        auto want = naive(dims, original, rhs);
        a.mul(b);
        assert(a.data == want && b.data == rhs && a.dim == dims);
        a.data = original;
        want = naive(dims, original, original);
        a.mul(a);
        assert(a.data == want);
        for(auto &x: a.data) {x = T::mod()-1;}
        b.data = a.data;
        want = naive(dims, a.data, b.data);
        a.mul(b);
        assert(a.data == want);
    }
    // Both mutable and const spans remain usable without an explicit template argument.
    big_vector<T> a{1, 2, 3, 4}, b{5, 6, 7, 8};
    auto x = subset_convolution(std::span(a), std::span(b));
    auto y = subset_convolution(std::span<T const>(a), std::span<T const>(b));
    assert(x == y && x == naive(big_vector<size_t>{2, 2}, a, b));
    // Exercise the rank cap and the largest supported ternary expansion with a closed-form oracle.
    std::vector<big_vector<size_t>> large{big_vector<size_t>(9, 2)};
    if(max_logn == 20) {
        large.push_back(big_vector<size_t>(20, 2));
        large.push_back({3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3});
    }
    for(auto const& dims: large) {
        fft::multivar<T> f(dims);
        std::ranges::fill(f.data, T::mod()-1);
        f.mul(f);
        for(size_t i = 0; i < f.N; i++) {
            size_t j = i;
            T want = 1;
            for(auto n: dims) {want *= T(j % n + 1); j /= n;}
            assert(f.data[i] == want);
        }
    }
}
int main() {
    check<modint<998244353>>();
    check<modint<1000000007>>();
    dynamic_modint<>::with_mod(998244353, [] {check<dynamic_modint<>>();});
    std::cout << "Multivariate products, squares, unit axes and const inputs passed under two primes and dynamic modint\n";
}
