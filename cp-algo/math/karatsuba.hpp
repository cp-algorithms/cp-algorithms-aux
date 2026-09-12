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
        auto a = &_a[0];
        auto b = &_b[0];
        auto c = &_c[0];
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
        using base = std::decay_t<decltype(a[0])>;
        uint64_t largest = base::mod() - 1;
        if(largest && largest > UINT64_MAX / N / largest) {
            base_conv<N>(a, b, c);
            return;
        }
        if constexpr (N % 4) {
            static_assert(N < 4);
            base_conv<N>(a, b, c);
            return;
        }
        alignas(32) uint64_t pr0[2 * N] = {}, pr1[2 * N] = {};
        alignas(32) uint64_t pr2[2 * N] = {}, pr3[2 * N] = {};
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
    // N is the input length and must be a power of two.
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

    namespace karatsuba_detail {
        template<class T>
        T constant(uint64_t x) {
            if constexpr(std::is_same_v<T, nimber::f2_64>) {
                T res{};
                res.r = x;
                return res;
            } else {
                return T(x);
            }
        }

        template<class T>
        T inverse(T a) {
            if constexpr(std::is_same_v<T, nimber::f2_64>) {
                return bpow(a, UINT64_MAX - 1, constant<T>(1));
            } else {
                // Euclid also supports composite moduli coprime to 2, 3, and 5.
                auto x = a.getr(), y = typename T::UInt(T::mod());
                typename T::Int2 u = 1, v = 0;
                while(y) {
                    auto q = x / y;
                    x -= q * y;
                    std::swap(x, y);
                    u -= q * v;
                    std::swap(u, v);
                }
                assert(x == 1);
                return T(u);
            }
        }

        template<class T>
        auto const& interpolation() {
            auto build = [] {
                std::array<std::array<T, 14>, 7> a{};
                // Values at 0, 1, 2, 3, 4, 5, and infinity.
                for(size_t i = 0; i < 6; i++) {
                    T x = constant<T>(i), p = constant<T>(1);
                    for(size_t j = 0; j < 7; j++) {a[i][j] = p; p *= x;}
                }
                a[6][6] = constant<T>(1);
                for(size_t i = 0; i < 7; i++) {a[i][i + 7] = constant<T>(1);}
                for(size_t i = 0; i < 7; i++) {
                    size_t pivot = i;
                    while(pivot < 7 && a[pivot][i] == T{}) {pivot++;}
                    assert(pivot < 7);
                    std::swap(a[i], a[pivot]);
                    auto inv = inverse(a[i][i]);
                    for(auto &x: a[i]) {x *= inv;}
                    for(size_t j = 0; j < 7; j++) {
                        if(j == i) {continue;}
                        auto q = a[j][i];
                        for(size_t k = 0; k < 14; k++) {a[j][k] -= q * a[i][k];}
                    }
                }
                std::array<std::array<T, 7>, 7> res{};
                for(size_t i = 0; i < 7; i++) {
                    for(size_t j = 0; j < 7; j++) {res[i][j] = a[i][j + 7];}
                }
                return res;
            };
            static thread_local auto matrix = build();
            if constexpr(modint_type<T>) {
                static thread_local auto modulus = T::mod();
                if(modulus != T::mod()) {matrix = build(); modulus = T::mod();}
            }
            return matrix;
        }

        // Toom-4 uses seven products for four blocks; Karatsuba handles the leaves.
        template<size_t N, class T>
        void mul(T const* a, T const* b, T* c) {
            if constexpr(N > (1 << 19)) {
                __builtin_unreachable(); // Same size bound as _karatsuba.
            } else if constexpr(N <= 4096) {
                std::fill_n(c, 2 * N, T{});
                _karatsuba<N>(a, b, c);
            } else {
                constexpr size_t h = N / 4;
                static big_vector<T> products(14 * h), av(h), bv(h);
                mul<h>(a, b, products.data());
                mul<h>(a + 3 * h, b + 3 * h, products.data() + 12 * h);
                for(size_t k = 1; k < 6; k++) {
                    auto point = constant<T>(k);
                    for(size_t i = 0; i < h; i++) {
                        av[i] = ((a[i + 3*h] * point + a[i + 2*h]) * point + a[i + h]) * point + a[i];
                        bv[i] = ((b[i + 3*h] * point + b[i + 2*h]) * point + b[i + h]) * point + b[i];
                    }
                    mul<h>(av.data(), bv.data(), products.data() + 2 * k * h);
                }
                auto const& matrix = interpolation<T>();
                std::fill_n(c, 2 * N, T{});
                for(size_t j = 0; j < 7; j++) {
                    for(size_t k = 0; k < 7; k++) {
                        auto x = matrix[j][k];
                        if(x == T{}) {continue;}
                        if(x == constant<T>(1)) {
                            for(size_t i = 0; i < 2*h - 1; i++) {c[j*h + i] += products[k*2*h + i];}
                        } else {
                            for(size_t i = 0; i < 2*h - 1; i++) {c[j*h + i] += x * products[k*2*h + i];}
                        }
                    }
                }
            }
        }
    }

    // Large field inputs use Toom-4 above the Karatsuba recursion.
    // Runtime wrapper that deduces N at compile time.
    // Resizes inputs to the next power of 2 and result to n + m - 1
    auto karatsuba(auto &a, auto &b) {
        using base = std::decay_t<decltype(a[0])>;
        auto n = std::size(a);
        auto m = std::size(b);
        if(!n || !m) {return big_vector<base>{};}
        auto N = std::bit_ceil(std::max(n, m));
        a.resize(N);
        b.resize(N);
        // The recursion reads the zero coefficient at index 2*N-1 when joining halves.
        big_vector<base> c(2 * N);
        with_bit_ceil(N, [&]<auto NN>() {
            if constexpr(std::is_same_v<base, nimber::f2_64>) {
                karatsuba_detail::mul<NN>(std::data(a), std::data(b), c.data());
                return;
            } else if constexpr(modint_type<base>) {
                if(base::mod() > 5 && base::mod() % 2 && base::mod() % 3 && base::mod() % 5) {
                    karatsuba_detail::mul<NN>(std::data(a), std::data(b), c.data());
                    return;
                }
            }
            _karatsuba<NN>(a, b, c);
        });
        c.resize(n + m - 1);
        return c;
    }
}

#endif // CP_ALGO_MATH_KARATSUBA_HPP
