#include "cp-algo/math/fft_simple.hpp"
#include <cassert>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <vector>

int main() {
    assert(std::fegetround() == FE_TONEAREST);
    auto check = [](double x) {
        // Constant real and imaginary parts have an exact FFT round trip.
        for(size_t k = 1; k <= 8; ++k) {
            std::vector<double> input(8, x);
            cp_algo::math::fft::dft_simple transform(input, 4);
            std::vector<int64_t> result(k);
            transform.recover_mod(result, k);
            for(auto value : result) {
                assert(value == int64_t(std::nearbyint(x)));
            }
        }
    };
    for(int i = -100; i <= 100; ++i) {
        double x = i * 0.5;
        check(x);
        check(std::nextafter(x, INFINITY));
        check(std::nextafter(x, -INFINITY));
    }
    for(double x : {0x1p51, 0x1p52, 0x1p62, 0x1p63}) {
        check(std::nextafter(x, 0.0));
        check(-std::nextafter(x, 0.0));
    }
    check(-0x1p63);
    // Compare signed convolution with exact integer arithmetic.
    for(int n : {1, 3, 4, 5, 17, 64, 65}) {
        std::vector<int64_t> a(n), b(n + 2), expected(2 * n + 1);
        for(int i = 0; i < n; ++i) a[i] = (i * 31 % 97) - 48;
        for(int j = 0; j < n + 2; ++j) b[j] = (j * 19 % 101) - 50;
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n + 2; ++j) expected[i + j] += a[i] * b[j];
        cp_algo::math::fft::conv_simple(a, b);
        assert(a == expected);
    }
}
