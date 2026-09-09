#ifndef CP_ALGO_LINALG_STRASSEN_HPP
#define CP_ALGO_LINALG_STRASSEN_HPP
#include "vector.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::linalg::impl {
    template<class row> constexpr bool use_strassen = false;
    template<auto mod> constexpr bool use_strassen<modint_vec<math::modint<mod>>> =
        mod > 1 && mod % 2 && mod < (1LL << 30);

    // Internal packed storage; the public matrix keeps its usual row representation.
    template<uint32_t mod>
    struct strassen_product {
        struct view {
            uint32_t *data;
            size_t stride;
            uint32_t* operator[](size_t i) const {return data + i * stride;}
        };
        static u64x4 mul(u64x4 a, u64x4 b) {
#ifdef __AVX2__
            return u64x4(_mm256_mul_epu32(__m256i(a), __m256i(b)));
#else
            return low32(a) * low32(b);
#endif
        }
        static u64x4 shrink(u64x4 x) {
            // x < 4*mod*2^32; reduce the high word modulo 2*mod.
            auto words = u32x8(x);
            auto bound = u32x8(u64x4() + (uint64_t(2) * mod << 32));
#ifdef __AVX2__
            return u64x4(_mm256_min_epu32(__m256i(words), __m256i(words - bound)));
#else
            return u64x4(words < words - bound ? words : words - bound);
#endif
        }
        template<size_t dim = 0>
        [[gnu::noinline]] static void leaf(view a, view b, view c, size_t n, size_t m, size_t k) {
            if constexpr(dim) {
                n = m = k = dim;
                a.stride = b.stride = c.stride = dim;
            }
            for(size_t i = 0; i < n; i += 4) {
                for(size_t j = 0; j < k; j += 8) {
                    u64x4 acc[4][2]{};
                    for(size_t first = 0; first < m; first += 8) {
#pragma GCC unroll 1
                        for(size_t z = first; z < first + 8; z++) {
                            u64x4 x;
                            std::memcpy(&x, b[z] + j, sizeof x);
                            u64x4 scales[4];
                            for(size_t t = 0; t < 4; t++) scales[t] = u64x4(u32x8() + a[i + t][z]);
                            for(size_t t = 0; t < 4; t++) acc[t][0] += mul(scales[t], x);
                            x >>= 32;
                            for(size_t t = 0; t < 4; t++) acc[t][1] += mul(scales[t], x);
                        }
                        for(auto &row: acc) for(auto &x: row) x = shrink(x);
                    }
                    for(size_t t = 0; t < 4; t++) {
                        constexpr uint32_t inv = math::inv2(uint32_t(-mod));
                        // x < 2*mod*2^32, so Montgomery reduction yields a value below 3*mod.
                        for(auto &x: acc[t]) x = montgomery_reduce(x, mod, inv);
                        u32x8 out = u32x8(acc[t][0] | (acc[t][1] << 32));
                        out = out < out - mod ? out : out - mod;
                        out = out < out - mod ? out : out - mod;
                        std::memcpy(c[i + t] + j, &out, sizeof out);
                    }
                }
            }
        }
        // c = a +/- b; c may alias a.
        template<bool subtract = false>
        static void combine(uint32_t *a, uint32_t *b, uint32_t *c, size_t n, size_t m) {
            for(size_t j = 0; j < n * m; j += 8) {
                u32x8 x, y;
                std::memcpy(&x, a + j, sizeof x);
                std::memcpy(&y, b + j, sizeof y);
                u32x8 z = subtract ? x + mod - y : x + y;
                z = z < z - mod ? z : z - mod;
                std::memcpy(c + j, &z, sizeof z);
            }
        }
        static bool can_split(size_t n, size_t m, size_t k) {
            return std::min({n, m, k}) > 64 && n % 16 == 0 && m % 16 == 0 && k % 16 == 0;
        }
        [[gnu::noinline]] static void multiply(uint32_t *a, uint32_t *b, uint32_t *c, size_t n, size_t m, size_t k,
                                              uint32_t *work) {
            // Every leaf dimension must remain a multiple of eight.
            if(!can_split(n, m, k)) {
                if(n == 64 && m == 64 && k == 64) leaf<64>({a, m}, {b, k}, {c, k}, n, m, k);
                else leaf({a, m}, {b, k}, {c, k}, n, m, k);
                return;
            }
            n /= 2; m /= 2; k /= 2;
            auto s = work, t = s + n * m, p = t + m * k;
            work += n * m + m * k + n * k;
            auto a00 = a, a01 = a + n * m, a10 = a + 2 * n * m, a11 = a + 3 * n * m;
            auto b00 = b, b01 = b + m * k, b10 = b + 2 * m * k, b11 = b + 3 * m * k;
            auto c00 = c, c01 = c + n * k, c10 = c + 2 * n * k, c11 = c + 3 * n * k;

            // Winograd's schedule uses seven products and fifteen additions.
            multiply(a00, b00, c11, n, m, k, work); // P1
            multiply(a01, b10, c00, n, m, k, work); // P2
            combine(c00, c11, c00, n, k); // C00 = P1 + P2

            combine(a10, a11, s, n, m); combine<true>(b01, b00, t, m, k); // S1, T1
            multiply(s, t, c01, n, m, k, work); // P5
            combine<true>(s, a00, s, n, m); combine<true>(b11, t, t, m, k); // S2, T2
            multiply(s, t, c10, n, m, k, work); // P6
            combine(c11, c10, c10, n, k); // U2 = P1 + P6

            combine<true>(a01, s, s, n, m); // S4
            multiply(s, b11, p, n, m, k, work); // P3
            combine(c10, c01, c11, n, k); // U4 = U2 + P5
            combine(c11, p, c01, n, k); // C01 = U4 + P3

            combine<true>(t, b10, t, m, k); // T4
            multiply(a11, t, p, n, m, k, work); // P4
            combine<true>(c10, p, c10, n, k);
            combine<true>(a00, a10, s, n, m); combine<true>(b11, b01, t, m, k); // S3, T3
            multiply(s, t, p, n, m, k, work); // P7
            combine(c10, p, c10, n, k); combine(c11, p, c11, n, k); // C10, C11
        }
        static u32x8 encode(u32x8 x) {
            // Scale by 2^32 with a fixed reciprocal; its quotient is off by at most one.
            constexpr uint32_t scale = (uint64_t(1) << 32) % mod;
            constexpr uint32_t quotient = (uint64_t(scale) << 32) / mod;
            auto packed = u64x4(x), q = u64x4() + quotient;
            auto lo = mul(packed, q) >> 32, hi = mul(packed >> 32, q);
            auto approx = u32x8(lo | (hi & (~uint64_t(0) << 32)));
            auto out = x * scale - approx * mod;
            return out < out - mod ? out : out - mod;
        }
        // Copy between matrix rows and contiguous recursive quadrants.
        template<bool unpack = false, class matrix>
        static void copy(matrix &a, uint32_t *ptr, size_t n, size_t m, size_t depth,
                         size_t row = 0, size_t col = 0) {
            if(depth) {
                n /= 2; m /= 2;
                for(size_t q = 0; q < 4; q++) {
                    copy<unpack>(a, ptr + q * n * m, n, m, depth - 1,
                                 row + q / 2 * n, col + q % 2 * m);
                }
                return;
            }
            if(col >= a.m()) return;
            size_t width = std::min(m, a.m() - col);
            for(size_t i = row; i < std::min(row + n, a.n()); i++) {
                auto data = a[i].data() + col;
                auto packed = ptr + (i - row) * m;
                for(size_t j = 0; j < width; j++) {
                    if constexpr(unpack) data[j].setr(packed[j]);
                    else packed[j] = uint32_t(data[j].getr());
                }
            }
        }
        template<class matrix>
        static matrix product(matrix const& a, matrix const& b) {
            auto pad = [](size_t x) {return (x + 31) / 32 * 32;};
            size_t n = pad(a.n()), m = pad(a.m()), k = pad(b.m());
            // Each recursion level needs a quarter as much scratch; children reuse it.
            size_t entries = n * m + m * k + n * k;
            // Small products avoid repeated mmap/madvise setup for temporary storage.
            std::vector<uint32_t> small;
            big_vector<uint32_t> large;
            size_t count = entries + entries / 3;
            auto ap = count < (1 << 21) ? (small.resize(count), small.data())
                                       : (large.resize(count), large.data());
            auto bp = ap + n * m, cp = bp + m * k;
            auto scratch = cp + n * k;
            // A, B, and C must use the same depth, including rectangular products.
            size_t depth = 0;
            for(size_t x = n, y = m, z = k; can_split(x, y, z); x /= 2, y /= 2, z /= 2) depth++;
            copy(a, ap, n, m, depth); copy(b, bp, m, k, depth);
            // Only A is scaled; each leaf's Montgomery reduction removes the factor.
            for(size_t i = 0; i < n * m; i += 8) {
                u32x8 x;
                std::memcpy(&x, ap + i, sizeof x);
                x = encode(x);
                std::memcpy(ap + i, &x, sizeof x);
            }
            multiply(ap, bp, cp, n, m, k, scratch);
            matrix res(a.n(), b.m());
            copy<true>(res, cp, n, k, depth);
            return res;
        }
    };
}
#pragma GCC pop_options
#endif // CP_ALGO_LINALG_STRASSEN_HPP
