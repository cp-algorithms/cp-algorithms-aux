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
            view at(size_t i, size_t j) const {return {data + i * stride + j, stride};}
        };
        static u64x4 mul(u64x4 a, u64x4 b) {
#ifdef __AVX2__
            return u64x4(_mm256_mul_epu32(__m256i(a), __m256i(b)));
#else
            return low32(a) * low32(b);
#endif
        }
        static u64x4 shrink(u64x4 x) {
            // x < 16*mod^2; subtraction puts the comparison in the signed range.
            auto y = x - uint64_t(8) * mod * mod;
            return i64x4(y) < 0 ? x : y;
        }
        [[gnu::noinline]] static void leaf(view a, view b, view c, size_t n, size_t m, size_t k) {
            for(size_t i = 0; i < n; i += 4) {
                for(size_t j = 0; j < k; j += 8) {
                    u64x4 acc[4][2]{};
                    for(size_t first = 0; first < m; first += 8) {
#pragma GCC unroll 1
                        for(size_t z = first; z < first + 8; z++) {
                            u64x4 x;
                            std::memcpy(&x, b[z] + j, sizeof x);
                            u64x4 y = x >> 32;
                            for(size_t t = 0; t < 4; t++) {
                                // Broadcast 32 bits: mul uses the low half of each 64-bit lane.
                                u64x4 scale = u64x4(u32x8() + a[i + t][z]);
                                acc[t][0] += mul(scale, x);
                                acc[t][1] += mul(scale, y);
                            }
                        }
                        for(auto &row: acc) for(auto &x: row) x = shrink(x);
                    }
                    for(size_t t = 0; t < 4; t++) {
                        constexpr uint32_t inv = math::inv2(uint32_t(-mod));
                        // x < 8*mod^2, so Montgomery reduction yields a value below 3*mod.
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
        static void combine(view a, view b, view c, size_t n, size_t m) {
            for(size_t i = 0; i < n; i++) for(size_t j = 0; j < m; j += 8) {
                u32x8 x, y;
                std::memcpy(&x, a[i] + j, sizeof x);
                std::memcpy(&y, b[i] + j, sizeof y);
                u32x8 z = subtract ? x + mod - y : x + y;
                z = z < z - mod ? z : z - mod;
                std::memcpy(c[i] + j, &z, sizeof z);
            }
        }
        [[gnu::noinline]] static void multiply(view a, view b, view c, size_t n, size_t m, size_t k) {
            // Every leaf dimension must remain a multiple of eight.
            if(std::min({n, m, k}) <= 64 || n % 16 || m % 16 || k % 16) {
                leaf(a, b, c, n, m, k);
                return;
            }
            n /= 2; m /= 2; k /= 2;
            big_vector<uint32_t> sb(n * m), tb(m * k), pb(n * k);
            view s{sb.data(), m}, t{tb.data(), k}, p{pb.data(), k};
            auto a00 = a, a01 = a.at(0, m), a10 = a.at(n, 0), a11 = a.at(n, m);
            auto b00 = b, b01 = b.at(0, k), b10 = b.at(m, 0), b11 = b.at(m, k);
            auto c00 = c, c01 = c.at(0, k), c10 = c.at(n, 0), c11 = c.at(n, k);
            auto copy = [&](view to) {for(size_t i = 0; i < n; i++) std::copy_n(p[i], k, to[i]);};

            combine(a00, a11, s, n, m); combine(b00, b11, t, m, k);
            multiply(s, t, p, n, m, k); // (A00 + A11)(B00 + B11)
            copy(c00); copy(c11);

            combine(a10, a11, s, n, m);
            multiply(s, b00, p, n, m, k); // (A10 + A11)B00
            copy(c10); combine<true>(c11, p, c11, n, k);

            combine<true>(b01, b11, t, m, k);
            multiply(a00, t, p, n, m, k); // A00(B01 - B11)
            copy(c01); combine(c11, p, c11, n, k);

            combine<true>(b10, b00, t, m, k);
            multiply(a11, t, p, n, m, k); // A11(B10 - B00)
            combine(c00, p, c00, n, k); combine(c10, p, c10, n, k);

            combine(a00, a01, s, n, m);
            multiply(s, b11, p, n, m, k); // (A00 + A01)B11
            combine<true>(c00, p, c00, n, k); combine(c01, p, c01, n, k);

            combine<true>(a10, a00, s, n, m); combine(b00, b01, t, m, k);
            multiply(s, t, p, n, m, k); // (A10 - A00)(B00 + B01)
            combine(c11, p, c11, n, k);

            combine<true>(a01, a11, s, n, m); combine(b10, b11, t, m, k);
            multiply(s, t, p, n, m, k); // (A01 - A11)(B10 + B11)
            combine(c00, p, c00, n, k);
        }
        template<class matrix>
        static matrix product(matrix const& a, matrix const& b) {
            auto pad = [](size_t x) {return (x + 31) / 32 * 32;};
            size_t n = pad(a.n()), m = pad(a.m()), k = pad(b.m());
            big_vector<uint32_t> ap(n * m), bp(m * k), cp(n * k);
            // Only A is scaled by 2^32; each leaf's Montgomery reduction removes it.
            for(size_t i = 0; i < a.n(); i++) for(size_t j = 0; j < a.m(); j++) {
                ap[i * m + j] = uint32_t((uint64_t(a[i][j].getr()) << 32) % mod);
            }
            for(size_t i = 0; i < b.n(); i++) for(size_t j = 0; j < b.m(); j++) {
                bp[i * k + j] = uint32_t(b[i][j].getr());
            }
            multiply({ap.data(), m}, {bp.data(), k}, {cp.data(), k}, n, m, k);
            matrix res(a.n(), b.m());
            for(size_t i = 0; i < res.n(); i++) for(size_t j = 0; j < res.m(); j++) {
                res[i][j].setr(cp[i * k + j]);
            }
            return res;
        }
    };
}
#pragma GCC pop_options
#endif // CP_ALGO_LINALG_STRASSEN_HPP
