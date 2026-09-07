#include "cp-algo/math/fft.hpp"
#include <random>
#include <iostream>
#include <utility>
using namespace cp_algo::math;
template<typename T> void check() {
    std::mt19937 rng(371);
    for(auto [small, large]: {std::pair{size_t(63), size_t(1 << 20)},
                             std::pair{size_t(64), size_t((1 << 20) - 1)},
                             std::pair{size_t(64), size_t(1 << 20)},
                             std::pair{size_t(1000), size_t((1 << 20) + 19)},
                             std::pair{size_t(4096), size_t((1 << 20) + 1)}}) {
        std::vector<T> a(large), b(small);
        for(auto &x: a) {x = rng() % T::mod();}
        for(auto &x: b) {x = rng() % T::mod();}
        size_t need = large + small - 1, length = std::bit_ceil(need);
        cp_algo::big_vector<T> expected(begin(a), end(a)), other(begin(b), end(b));
        expected.resize(length); other.resize(length);
        // Independent existing cyclic backend, bypassing unbalanced dispatch.
        fft::cyclic_mul(expected, other, length, large <= length / 2);
        expected.resize(need);
        auto original = a;
        fft::mul(a, std::as_const(b));
        assert(std::ranges::equal(a, expected));
        if(small == 1000) {
            fft::mul(b, std::as_const(original));
            assert(std::ranges::equal(b, expected));
        }
    }
    size_t n = (1 << 20) + 17, m = 100;
    cp_algo::big_vector<T> a(n, 1), b(m, 1);
    fft::mul(a, b);
    for(size_t i = 0; i < a.size(); i++) {
        assert(a[i] == T(std::min({i + 1, n, m, n + m - 1 - i})));
    }
}
int main() {
    check<modint<998244353>>();
    check<modint<1000000007>>();
    std::cout << "Unbalanced convolution properties passed under two primes\n";
}
