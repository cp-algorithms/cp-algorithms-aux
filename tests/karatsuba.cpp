#include "cp-algo/math/karatsuba.hpp"
#include <random>

using namespace cp_algo;
using namespace cp_algo::math;

template<class T>
T coefficient(uint64_t x) {
    if constexpr(std::is_same_v<T, nimber::f2_64>) {
        T res{};
        res.r = x;
        return res;
    } else if constexpr(modint_type<T>) {
        return T(x % T::mod());
    } else {
        return T(x % 7);
    }
}

template<class T>
void check(size_t n, size_t m, bool square = false) {
    std::mt19937_64 rng(100 * n + m);
    std::vector<T> a(n), b(m);
    for(auto &x: a) {x = coefficient<T>(rng());}
    for(auto &x: b) {x = coefficient<T>(rng());}
    if(square) {b = a;}
    big_vector<T> expected(n && m ? n + m - 1 : 0);
    if(std::min(n, m) < 80) {
        for(size_t i = 0; i < n; i++) {
            for(size_t j = 0; j < m; j++) {expected[i + j] += a[i] * b[j];}
        }
    } else {
        // Compare Toom interpolation against the original three-product recursion.
        size_t N = std::bit_ceil(std::max(n, m));
        auto x = a, y = b;
        x.resize(N); y.resize(N);
        expected.resize(2 * N);
        with_bit_ceil(N, [&]<auto NN>() {_karatsuba<NN>(x, y, expected);});
        expected.resize(n + m - 1);
    }
    auto actual = square ? karatsuba(a, a) : karatsuba(a, b);
    assert(actual == expected);
}

template<class T>
void boundaries() {
    for(size_t n: {0, 1, 2, 3, 7, 8, 9, 15, 16, 17, 33, 65}) {
        for(size_t m: {0, 1, 3, 8, 17, 65}) {check<T>(n, m);}
    }
    for(size_t n: {4095, 4096, 4097, 8192, 8193}) {
        check<T>(n, n - 13);
        check<T>(n, 17);
    }
    check<T>(8193, 8193, true);
}

int main() {
    boundaries<modint<1000000007>>();
    boundaries<nimber::f2_64>();
    boundaries<modint<2147483647>>(); // Leaf sums must not overflow uint64_t.
    check<modint<1>>(8193, 8193);
    check<modint<2>>(8193, 8193);
    check<modint<15>>(8193, 8193); // Noninvertible interpolation denominators.
    check<modint<49>>(8193, 8193); // Composite, but all denominators invertible.
    check<int64_t>(65, 65);
    for(int modulus: {1000000007, 49, 37, 15, 1000000007}) {
        dynamic_modint<>::with_mod(modulus, [] {check<dynamic_modint<>>(8193, 8179);});
    }
    check<modint<2305843009213693951LL>>(4097, 4097);
    std::cout << "Karatsuba and Toom-4 properties passed\n";
}
