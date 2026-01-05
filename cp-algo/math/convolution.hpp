#ifndef CP_ALGO_MATH_CONVOLUTION_HPP
#define CP_ALGO_MATH_CONVOLUTION_HPP
#include "fft.hpp"
#include "cvector.hpp"
#include <vector>
#include <algorithm>
#include <bit>
#include <type_traits>
#include <ranges>

namespace cp_algo::math {

// Convolution limited to the first `need` coefficients.
// Writes the result into `a`; performs in-place when possible (modint path).
template<class VecA, class VecB>
void convolution_prefix(VecA& a, VecB const& b, size_t need) {
    using T = std::decay_t<decltype(a[0])>;
    if constexpr (modint_type<T>) {
        // Use NTT-based truncated multiplication. Works in-place on `a`.
        fft::mul_truncate(a, b, need);
        return;
    }
    size_t na = std::min(need, std::size(a));
    size_t nb = std::min(need, std::size(b));
    a.resize(na);
    auto bv = b | std::views::take(nb);

    if(na == 0 || nb == 0) {
        a.clear();
        return;
    }
    if constexpr (std::is_same_v<T, fft::point>) {
        size_t conv_len = na + nb - 1;
        size_t n = std::bit_ceil(conv_len);
        n = std::max(n, (size_t)fft::flen);
        fft::cvector A(n), B(n);
        for(size_t i = 0; i < na; i++) {
            A.set(i, a[i]);
        }
        for(size_t i = 0; i < nb; i++) {
            B.set(i, bv[i]);
        }
        A.fft();
        B.fft();
        A.dot(B);
        A.ifft();
        a.assign(need, T(0));
        for(size_t i = 0; i < std::min(need, conv_len); i++) {
            a[i] = A.template get<fft::point>(i);
        }
    } else if constexpr (std::is_same_v<T, fft::ftype>) {
        // Imaginary-cyclic convolution modulo x^n-i to compute acyclic convolution
        // Represents polynomials as point(a[i], a[i+n]) to work in x^n-i basis
        size_t conv_len = na + nb - 1;
        size_t n = std::bit_ceil(conv_len) / 2;
        n = std::max(n, (size_t)fft::flen);
        
        fft::cvector A(n), B(n);
        // Pack as modulo x^n-i: A[i] = point(a[i], a[i+n])
        for(size_t i = 0; i < std::min(n, na); i++) {
            fft::ftype re = a[i], im = 0;
            if(i + n < na) im = a[i + n];
            A.set(i, fft::point(re, im));
        }
        for(size_t i = 0; i < std::min(n, nb); i++) {
            fft::ftype re = bv[i], im = 0;
            if(i + n < nb) im = bv[i + n];
            B.set(i, fft::point(re, im));
        }
        A.fft();
        B.fft();
        A.dot(B);
        A.ifft();
        a.assign(2 * n, T(0));
        for(size_t i = 0; i < n; i++) {
            auto v = A.template get<fft::point>(i);
            a[i] = v.real();
            a[i + n] = v.imag();
        }
        a.resize(need);
    } else {
        // Generic fallback: use simple O(n^2) convolution from fft utilities.
        fft::mul_slow(a, bv, need);
    }
}

} // namespace cp_algo::math

#endif // CP_ALGO_MATH_CONVOLUTION_HPP
