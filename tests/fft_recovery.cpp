#include "cp-algo/math/fft_simple.hpp"
#include <cassert>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <vector>

int main() {
    int original_mode = std::fegetround();
    for(int mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
        assert(std::fesetround(mode) == 0);
        auto check = [](double x) {
            // Constant real and imaginary parts have an exact FFT round trip.
            std::vector<double> input(8, x);
            cp_algo::math::fft::dft_simple transform(input, 4);
            std::vector<int64_t> result(8);
            transform.recover_mod(result, result.size());
            for(auto value : result) {
                assert(value == std::llround(x));
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
    }
    assert(std::fesetround(original_mode) == 0);
}
