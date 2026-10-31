#ifndef CP_ALGO_MATH_FFT_SIMPLE_HPP
#define CP_ALGO_MATH_FFT_SIMPLE_HPP
#include "../random/rng.hpp"
#include "../math/common.hpp"
#include "../math/cvector.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::fft {
    struct dft_simple {
        cp_algo::math::fft::cvector cv;

        dft_simple(auto const& a, size_t n): cv(n) {
            for(size_t i = 0; i < std::min(std::size(a), n); i++) {
                real(cv.at(i))[i % 4] = ftype(a[i]);
                imag(cv.at(i))[i % 4] = ftype(i + n < std::size(a) ? a[i + n] : 0);
            }
            checkpoint("dft64 init");
            cv.fft();
        }

        void dot(dft_simple const& t) {
            cv.dot(t.cv);
        }

        void recover_mod(auto &res, size_t k) {
            cv.ifft();
            size_t n = cv.size();
            // Finite coefficients whose nearest integer fits int64_t.
            // Truncation and the exact fractional remainder preserve llround's
            // ties-away-from-zero rule independently of the rounding mode.
            auto recover_integer = [](double x) -> int64_t {
                int64_t q = static_cast<int64_t>(x);
                double f = x - static_cast<double>(q);
                return q + (f >= 0.5) - (f <= -0.5);
            };
            for(size_t i = 0; i < std::min(k, n); i++) {
                res[i] = recover_integer(real(cv.get(i)));
            }
            for(size_t i = n; i < k; i++) {
                res[i] = recover_integer(imag(cv.get(i - n)));
            }
            cp_algo::checkpoint("recover mod");
        }
    };

    // Multiplies a and b, assuming perfect precision and no overflow
    void conv_simple(auto& a, auto const& b) {
        if (empty(a) || empty(b)) {
            a.clear();
            return;
        }
        size_t n = a.size(), m = b.size();
        size_t N = std::max(flen, std::bit_ceil(n + m - 1) / 2);
        dft_simple A(a, N), B(b, N);
        A.dot(B);
        a.resize(n + m - 1);
        A.recover_mod(a, n + m - 1);
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_FFT_SIMPLE_HPP
