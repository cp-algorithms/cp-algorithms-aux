#include "cp-algo/math/fft64.hpp"
#include <cassert>
#include <cfenv>
#include <cstdint>
#include <random>
#include <vector>

int main() {
    assert(std::fegetround() == FE_TONEAREST);
    std::mt19937_64 rng(42);
    for(uint64_t factor : {uint64_t(1), uint64_t(3), UINT64_MAX, rng() | 1}) {
        cp_algo::math::fft::dft64::_init = true;
        cp_algo::math::fft::dft64::factor = factor;
        cp_algo::math::fft::dft64::ifactor = cp_algo::math::inv2(factor);
        for(int n : {1, 3, 4, 5, 17, 64, 65, 129}) {
            for(int m : {1, 4, 9, 66}) {
                for(bool dense : {false, true}) {
                    std::vector<uint64_t> a(n), b(m), expected(n + m - 1);
                    for(auto &x : a) x = dense ? UINT64_MAX : rng();
                    for(auto &x : b) x = dense ? UINT64_MAX : rng();
                    for(int i = 0; i < n; ++i)
                        for(int j = 0; j < m; ++j) expected[i + j] += a[i] * b[j];
                    cp_algo::math::fft::conv64(a, b);
                    assert(a == expected);
                }
            }
        }
    }
}
