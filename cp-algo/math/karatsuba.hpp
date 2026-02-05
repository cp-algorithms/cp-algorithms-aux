#ifndef CP_ALGO_MATH_KARATSUBA_HPP
#define CP_ALGO_MATH_KARATSUBA_HPP
#include "../number_theory/nimber.hpp"
#include "../number_theory/modint.hpp"
#include "../util/big_alloc.hpp"
#include "../util/bit.hpp"
#include <vector>
#include <bit>
#include <cstdint>
#include <span>

namespace cp_algo::math {
    constexpr size_t NN = 8;

    template<auto N>
    void base_conv(auto &&_a, auto &&_b, auto &&_c) {
        auto a = std::assume_aligned<32>(&_a[0]);
        auto b = std::assume_aligned<32>(&_b[0]);
        auto c = std::assume_aligned<32>(&_c[0]);
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                c[i + j] += a[i] * b[j];
            }
        }
    }

    // Optimized base case for F2_64: uses 256-bit VPCLMULQDQ
    // Computes 4 products per iteration
    template<size_t N>
    [[gnu::target("avx2,vpclmulqdq")]]
    void base_conv_f2_64(auto &&a, auto &&b, auto &&c) {
        if constexpr (N % 2) {
            static_assert(N < 2);
            base_conv<N>(a, b, c);
            return;
        }
        alignas(32) __m128i pr0[2 * N] = {};
        alignas(32) __m128i pr1[2 * N] = {};
        
        for (size_t i = 0; i + 1 < N; i += 2) {
            auto va = (__m256i)u64x4{a[i], 0, a[i + 1], 0};
            for (size_t j = 0; j + 1 < N; j += 2) {
                auto vb = (__m256i)u64x4{b[j], b[j + 1], b[j], b[j + 1]};
                (__m256i&)pr0[i + j] ^= _mm256_clmulepi64_epi128(va, vb, 0);
                (__m256i&)pr1[i + j] ^= _mm256_clmulepi64_epi128(va, vb, 16);
            }
        }
        c[0].r = nimber::reduce_mod(pr0[0]);
        for (size_t i = 1; i < 2 * N - 1; i++) {
            c[i].r ^= nimber::reduce_mod(pr0[i] ^ pr1[i - 1]);
        }
    }

    template<auto N>
    void base_conv_modint(auto &&a, auto &&b, auto &&c) {
        if constexpr (N % 4) {
            static_assert(N < 4);
            base_conv<N>(a, b, c);
            return;
        }
        alignas(32) uint64_t pr0[2 * N] = {}, pr1[2 * N] = {};
        alignas(32) uint64_t pr2[2 * N] = {}, pr3[2 * N] = {};
        using base = std::decay_t<decltype(a[0])>;
        for (size_t i = 0; i < N; i += 4) {
            auto va0 = __m256i() + a[i].getr();
            auto va1 = __m256i() + a[i + 1].getr();
            auto va2 = __m256i() + a[i + 2].getr();
            auto va3 = __m256i() + a[i + 3].getr();
            size_t j = 0;
            for (; j + 3 < N; j += 4) {
                auto vb = (__m256i)u64x4{
                    b[j].getr(), b[j + 1].getr(), b[j + 2].getr(), b[j + 3].getr()
                };
                (__m256i&)pr0[i + j] += _mm256_mul_epu32(va0, vb);
                (__m256i&)pr1[i + j] += _mm256_mul_epu32(va1, vb);
                (__m256i&)pr2[i + j] += _mm256_mul_epu32(va2, vb);
                (__m256i&)pr3[i + j] += _mm256_mul_epu32(va3, vb);
            }
        }
        for (size_t i = 0; i < 2 * N - 1; i++) {
            if (i > 0) {
                pr2[i] += pr3[i - 1];
                pr1[i] += pr2[i - 1];
                pr0[i] += pr1[i - 1];
            }
            c[i].setr((typename base::UInt)(pr0[i] % base::mod()));
        }
    }

    // Generic Karatsuba multiplication algorithm for polynomials
    // Template parameters:
    //   N - Size of input arrays (must be power of 2)
    //   Add, Sub, Mul - Operations for addition, subtraction, and coefficient multiplication
    template<auto N>
    void _karatsuba(auto &&a, auto &&b, auto &&c) {
        [[gnu::assume(N <= 1<<19)]];
        using base = std::decay_t<decltype(a[0])>;
        if constexpr (N <= NN) {
            if constexpr (std::is_same_v<base, nimber::f2_64>) {
                base_conv_f2_64<N>(a, b, c);
            } else if constexpr (modint_type<base>) {
                base_conv_modint<N>(a, b, c);
            } else {
                base_conv<N>(a, b, c);
            }
        } else {
            constexpr auto h = N / 2;
            auto a0 = &a[0], a1 = a0 + h, b0 = &b[0], b1 = b0 + h;
            auto c0 = &c[0], c1 = c0 + h, c2 = c0 + 2 * h;
            _karatsuba<h>(a0, b0, c0);
            _karatsuba<h>(a1, b1, c2);
            static big_vector<base> buf(4 * h);
            auto f = &buf[0];
            auto sum_a = f + 2 * h, sum_b = f + 3 * h;
            for (size_t i = 0; i < h; i++) {
                sum_a[i] = a0[i] + a1[i];
                sum_b[i] = b0[i] + b1[i];
            }
            memset(f, 0, sizeof(base) * 2 * h);
            _karatsuba<h>(sum_a, sum_b, f);
            for(size_t i = 0; i < h; i++) {
                auto A = c0[i], &B = c1[i], &C = c2[i], D = c2[i + h];
                auto BC = B - C;
                B = BC + f[i] - A;
                C = f[i + h] - D - BC;
            }
        }
    }

    // Runtime wrapper that deduces N at compile time
    // Resizes inputs to the next power of 2 and result to n + m - 1
    auto karatsuba(auto &a, auto &b) {
        auto n = std::size(a);
        auto m = std::size(b);
        auto N = std::bit_ceil(std::max(n, m));
        a.resize(N);
        b.resize(N);
        using base = std::decay_t<decltype(a[0])>;
        big_vector<base> c(2 * N - 1);
        with_bit_ceil(N, [&]<auto NN>() {
            _karatsuba<NN>(a, b, c);
        });
        c.resize(n + m - 1);
        return c;
    }
}

#endif // CP_ALGO_MATH_KARATSUBA_HPP
