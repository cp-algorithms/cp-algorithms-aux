#include "cp-algo/math/fft.hpp"
#include <random>
#include <iostream>
using namespace cp_algo;
using namespace cp_algo::math;

template<int mod> void square_aliases() {
    using T = modint<mod>;
    std::mt19937 rng(42);
    for(size_t n: {1, 31, 63, 64, 65, 127, 129, 257, 1025, 65536, 65537,
                  524288, 524289, 1048576, 1048583}) {
        big_vector<T> a(n);
        for(size_t i = 0; i < n; i++) {a[i] = i % 2 ? T(-1) : T(1);}
        for(bool immutable: {false, true}) {
            auto b = a;
            if(immutable) {fft::mul(b, std::as_const(b));}
            else {fft::mul(b, b);}
            assert(b.size() == 2 * n - 1);
            // Exact linear-time reference on both sides of the large cutoff.
            for(size_t i = 0; i < b.size(); i++) {
                T want = i < n ? i + 1 : 2 * n - 1 - i;
                if(i % 2) {want = -want;}
                assert(b[i] == want);
            }
        }
        if(n <= 1025 || n == 65537 || n == 524289 || n == 1048583) {
            for(auto &x: a) {x = rng() % mod;}
            // Independent operands use the general multiplication path.
            auto expected = a, rhs = a;
            fft::mul(expected, rhs);
            auto actual = a;
            fft::mul(actual, actual);
            assert(actual == expected);
            if(n <= 1025 || n == 524289) {
                for(size_t k: {n, 2 * n - 1, 2 * n + 7}) {
                    actual = a;
                    auto want = expected;
                    want.resize(k);
                    fft::mul_truncate(actual, actual, k);
                    assert(actual == want);
                }
            }
        }
    }
}
int main() {
    cp_algo::random::gen.seed(42);
    square_aliases<998244353>();
    square_aliases<1000000007>();
    square_aliases<17>();
    square_aliases<65537>();
    std::cout << "Self-squares, aliased truncation and exact large references passed under four primes\n";
}
