#ifndef CP_ALGO_MATH_RING_HPP
#define CP_ALGO_MATH_RING_HPP
#include "../number_theory/discrete_sqrt.hpp"
#include "../number_theory/primality.hpp"
#include "../number_theory/modint.hpp"
#include "../util/checkpoint.hpp"
#include "../random/rng.hpp"
#include "cvector.hpp"
#include <cstring>
#include <ranges>
#include <tuple>
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::fft {
    size_t com_size(size_t as, size_t bs) {
        if(!as || !bs) {
            return 0;
        }
        return std::max(flen, std::bit_ceil(as + bs - 1) / 2);
    }
    // Whether a product of these sizes exceeds half of its padded length by a tail short enough
    // to be corrected naively, which lets mul_truncate halve its transform.
    constexpr size_t short_tail = 32;
    bool has_short_tail(size_t as, size_t bs) {
        size_t n = com_size(as, bs);
        return as + bs - 1 - n <= short_tail && as <= n && bs <= n;
    }
    // Convolution through an imaginary quadratic ring, for a prime modulus p < 2^31.
    // Let d be the smallest positive integer such that -d is a quadratic residue (d = 1 exactly
    // when p = 1 mod 4) and root^2 = -d. A residue x is written as re + im*sqrt(-d) with
    // re + im*root = x, (re, im) reduced by the lattice of representations of zero under the norm
    // re^2 + d*im^2, so that |re| and sqrt(d)*|im| are of order sqrt(p), and is embedded as the
    // complex number re + i*sqrt(d)*im. The product is then one complex convolution instead of
    // three real ones. For d = 1 (Gaussian integers) it is computed modulo x^n - i and x^n + i
    // (the latter through conjugation) and recombined modulo p, which needs i to lie in the ring;
    // for any other d it is a single transform of the full product length.
    constexpr uint64_t pow_mod(uint64_t x, uint64_t e, uint64_t p) {
        uint64_t res = 1;
        for(x %= p; e; e >>= 1, x = x * x % p) {if(e & 1) {res = res * x % p;}}
        return res;
    }
    // Smallest d > 0 with -d a quadratic residue modulo the odd prime p < 2^31, or 0.
    constexpr uint32_t least_negated_residue(uint64_t p) {
        for(uint32_t d = 1; d < 256 && d < p; d++) {
            if(pow_mod(p - d, (p - 1) / 2, p) == 1) {return d;}
        }
        return 0;
    }
    template<modint_type base>
    struct quadratic {
        static inline bool ready = false, available = false;
        static inline uint32_t d = 0, root = 0, prime = 0;
        // (a, b) and (a2, b2): reduced basis of the pairs (re, im) that represent zero.
        static inline int32_t a = 0, b = 0, a2 = 0, b2 = 0;
        // A fixed modulus settles d at compile time, so only the path in use is instantiated.
        static constexpr bool fixed_mod = requires {typename std::bool_constant<(base::mod(), true)>;};
        static constexpr uint32_t fixed_d = [] {
            if constexpr(fixed_mod) {return base::mod() > 2 && base::mod() % 2 ? least_negated_residue(base::mod()) : 0u;}
            else {return 0u;}
        }();
        static void init() {
            if(ready && prime == uint32_t(base::mod())) {return;}
            ready = true; available = false;
            uint64_t p = base::mod();
            prime = uint32_t(p);
            if(p < 3 || p >= (uint64_t(1) << 31) || !is_prime(p)) {return;}
            d = least_negated_residue(p);
            auto s = d ? cp_algo::math::sqrt(base(-int64_t(d))) : std::nullopt;
            if(!s) {return;}
            root = s->getr();
            // Lagrange reduction of (p, 0), (-root, 1) under the norm x^2 + d*y^2.
            using wide = __int128;
            wide x1 = p, y1 = 0, x2 = -wide(root), y2 = 1;
            auto norm = [&](wide x, wide y) {return x * x + d * y * y;};
            while(true) {
                if(norm(x1, y1) > norm(x2, y2)) {std::swap(x1, x2); std::swap(y1, y2);}
                wide dot = x1 * x2 + d * y1 * y2, len = norm(x1, y1);
                wide m = (2 * dot + (dot >= 0 ? len : -len)) / (2 * len);
                if(m == 0 || norm(x2 - m * x1, y2 - m * y1) >= norm(x2, y2)) {break;}
                x2 -= m * x1; y2 -= m * y1;
            }
            a = int32_t(x1); b = int32_t(y1);
            // In the Gaussian integers i*(a, b) = (-b, a) lies in the lattice as well.
            if(d == 1) {a2 = b; b2 = -a;}
            else {a2 = int32_t(x2); b2 = int32_t(y2);}
            assert(base(a) + base(b) * base(root) == base(0) && base(a2) + base(b2) * base(root) == base(0));
            assert(std::abs(int64_t(a) * b2 - int64_t(a2) * b) == int64_t(p));
            available = true;
        }
        // Points per transform: two branches of half the padded product length for d = 1, one
        // transform of the whole length otherwise. A short tail beyond that is computed naively
        // and taken back out, as in mul_truncate, which halves the transform: over the residues
        // for d = 1, where the branches together work modulo x^(2n) + 1, and in ring coordinates
        // for the single transform, which works modulo x^n - i.
        static size_t length(size_t as, size_t bs) {
            return (d == 1 ? com_size(as, bs) : 2 * com_size(as, bs)) >> has_short_tail(as, bs);
        }
        // Whether this path is used for operands of these sizes. With n = com_size(as, bs) it
        // takes six transforms of n points (three of 2n for d > 1), four for a square, where the
        // split representation takes seven and five, so it is preferred whenever it is exact,
        // unless a transform exceeds 2^24 points, the only length with a kernel tiled for data
        // outside the caches (which in turn does not let an operand wrap around).
        // Exactness: the largest rounding error measured is about
        // 0.003 * sqrt(need * d / 2^22) * p / 2^30, so the bound below keeps it under 1/16.
        static bool usable(size_t as, size_t bs) {
            init();
            if(!available || std::min(as, bs) < size_t(magic)) {return false;}
            size_t need = as + bs - 1, n = length(as, bs);
            if(n > (1 << 24) || (n == (1 << 24) && d == 1 && std::max(as, bs) > n)) {return false;}
            using wide = unsigned __int128;
            return wide(need) * d * base::mod() * base::mod() <= wide(1) << 90;
        }
        // Lattice constants as doubles, hoisted out of the hot loops (the statics are 32-bit
        // integers and could alias the 32-bit output stream).
        struct lattice {
            // (c, e): the basis vector with the larger second coordinate, used to shrink im.
            double a, b, a2, b2, to_q, to_t, c, e, inv_e, root, scale, inv_scale, p;
            lattice(): a(quadratic::a), b(quadratic::b), a2(quadratic::a2), b2(quadratic::b2), p(base::mod()) {
                double det = a * b2 - a2 * b;
                to_q = b2 / det; to_t = -b / det;
                bool first = std::abs(b) >= std::abs(b2);
                c = first ? a : a2; e = first ? b : b2; inv_e = 1.0 / e;
                root = quadratic::root;
                scale = std::sqrt(double(d)); inv_scale = 1.0 / scale;
            }
        };
        // Map a rounded coordinate pair back to the residue re + root*im.
        template<bool unit>
        static u32x4 project(vpoint value, bool negative, lattice const& L) {
            const double p = L.p;
            auto R = round(real(value)), I = negative ? -imag(value) : imag(value);
            if constexpr(unit) {I = round(I);}
            else {I = round(I * L.inv_scale);}
            auto q = round(I * L.inv_e);
            auto U = vftype(_mm256_fnmadd_pd(__m256d(q), _mm256_set1_pd(L.c), __m256d(R)));
            auto V = vftype(_mm256_fnmadd_pd(__m256d(q), _mm256_set1_pd(L.e), __m256d(I)));
            auto h = vftype(_mm256_fmadd_pd(__m256d(V), _mm256_set1_pd(L.root), __m256d(U)));
            q = round(h * (1.0 / p));
            auto out = vftype(_mm256_fnmadd_pd(__m256d(q), _mm256_set1_pd(p), __m256d(h)));
            out = out < 0 ? out + p : out;
            return u32x4(_mm256_cvttpd_epi32(__m256d(out)));
        }
        // Lift residues to ring coordinates with stochastic rounding.
        // The xorshift stream is advanced once per 8 residues and restarted from the same
        // seed for both branches, so both see the same representatives.
        // With stream set the spectrum is written with non-temporal stores and left for
        // cache_product to transform: for n = 2^24 it is far larger than the caches and is next
        // read by a separate pass, so this saves the read-for-ownership of every destination line.
        // Coefficients from n on (d = 1 only) wrap around: x^n = i in the branch modulo x^n - i,
        // and the conjugated branch modulo x^n + i sees conj(-i * z) = i * conj(z).
        template<bool stream, bool unit>
        static void fill(cvector& c, auto const& x, auto const& upper, size_t n, bool negative, u64x4 state, bool transform = true) {
            const lattice L;
            auto const* src = reinterpret_cast<const uint32_t*>(std::data(x));
            size_t count = std::size(x);
            assert(count <= n && (std::empty(upper) || (unit && !stream && count == n)));
            c.r.resize(n / flen);
            auto* dst = reinterpret_cast<double*>(c.r.data());
            u32x8 words{};
            auto emit = [&]<bool wrap = false>(size_t i, i32x4 bits) __attribute__((always_inline)) {
                auto v = __builtin_convertvector(bits, vftype);
                if(i % 8 == 0) {state ^= state << 13; state ^= state >> 7; state ^= state << 17; words = u32x8(state);}
                i32x4 small = i % 8 ? i32x4(__builtin_shufflevector(words, words, 1, 3, 5, 7)) : i32x4(__builtin_shufflevector(words, words, 0, 2, 4, 6));
                auto noise = __builtin_convertvector(small, vftype) * 0x1p-32;
                auto q = round(v * L.to_q + noise), t = round(v * L.to_t + noise);
                auto re = v - q * L.a - t * L.a2, im = -(q * L.b) - t * L.b2;
                if constexpr(!unit) {im = im * L.scale;}
                if constexpr(stream) {
                    _mm256_stream_pd(dst + 2 * i, __m256d(re));
                    _mm256_stream_pd(dst + 2 * i + flen, __m256d(negative ? -im : im));
                } else if constexpr(wrap) {
                    c.r[(i - n) / flen] += vpoint{negative ? im : -im, re};
                } else {
                    c.r[i / flen] = vpoint{re, negative ? -im : im};
                }
            };
            size_t full = count / flen * flen;
            for(size_t i = 0; i < full; i += flen) {
                i32x4 bits;
                std::memcpy(&bits, src + i, sizeof(bits));
                emit(i, bits);
            }
            if(full < count) {
                i32x4 bits{};
                for(size_t j = full; j < count; j++) {bits[j - full] = int32_t(src[j]);}
                emit(full, bits);
            }
            if constexpr(stream) {_mm_sfence();}
            std::fill(c.r.begin() + (count + flen - 1) / flen, c.r.end(), vpoint{});
            auto const* wrapped = reinterpret_cast<const uint32_t*>(std::data(upper));
            for(size_t i = 0; i < std::size(upper); i += flen) {
                i32x4 bits{};
                for(size_t j = i; j < std::min(i + flen, std::size(upper)); j++) {bits[j - i] = int32_t(wrapped[j]);}
                emit.template operator()<true>(n + i, bits);
            }
            checkpoint("quadratic init");
            if constexpr(!stream) {if(transform) {c.forward();}}
        }
        // Lift a range of residues to ring coordinates, one point per coefficient, for the
        // reusable transforms below. Unlike fill it reads through the modint interface, so it
        // serves views and Montgomery storage alike, and it neither wraps nor streams.
        template<bool unit>
        static void lift(cvector& c, auto const& x, size_t n, bool negative, u64x4 state) {
            const lattice L;
            size_t total = std::size(x), count = std::min(n, total);
            assert(total <= 2 * n && (d == 1 || total <= n));
            c.r.resize(n / flen);
            u32x8 words{};
            size_t step = 0;
            auto emit = [&]<bool wrap>(size_t i, i32x4 bits) __attribute__((always_inline)) {
                auto v = __builtin_convertvector(bits, vftype);
                if(step % 2 == 0) {state ^= state << 13; state ^= state >> 7; state ^= state << 17; words = u32x8(state);}
                i32x4 small = step++ % 2 ? i32x4(__builtin_shufflevector(words, words, 1, 3, 5, 7))
                                         : i32x4(__builtin_shufflevector(words, words, 0, 2, 4, 6));
                auto noise = __builtin_convertvector(small, vftype) * 0x1p-32;
                auto q = round(v * L.to_q + noise), t = round(v * L.to_t + noise);
                auto re = v - q * L.a - t * L.a2, im = -(q * L.b) - t * L.b2;
                if constexpr(!unit) {im = im * L.scale;}
                // Coefficients from n on wrap around: x^n is i in this branch, -i in the other.
                if constexpr(wrap) {c.r[(i - n) / flen] += vpoint{negative ? im : -im, re};}
                else {c.r[i / flen] = vpoint{re, negative ? -im : im};}
            };
            // Residues are raw words exactly when the modulus is a compile-time constant; wider
            // storage than the residue needs is narrowed on the way in.
            constexpr bool plain = fixed_mod && std::ranges::contiguous_range<std::decay_t<decltype(x)>>;
            constexpr bool raw = plain && sizeof(base) == 4;
            constexpr bool wide = plain && sizeof(base) == 8;
            auto load = [&](size_t i, size_t upto) {
                i32x4 bits{};
                if constexpr(raw || wide) {
                    if(i + flen <= upto) {
                        if constexpr(raw) {
                            std::memcpy(&bits, reinterpret_cast<const uint32_t*>(std::data(x)) + i, sizeof(bits));
                        } else {
                            u64x4 words;
                            std::memcpy(&words, reinterpret_cast<const uint64_t*>(std::data(x)) + i, sizeof(words));
                            bits = __builtin_convertvector(words, i32x4);
                        }
                        return bits;
                    }
                }
                for(size_t j = i; j < std::min(i + flen, upto); j++) {bits[j - i] = int32_t(x[j].getr());}
                return bits;
            };
            for(size_t i = 0; i < count; i += flen) {emit.template operator()<false>(i, load(i, count));}
            std::fill(c.r.begin() + (count + flen - 1) / flen, c.r.end(), vpoint{});
            for(size_t i = n; i < total; i += flen) {emit.template operator()<true>(i, load(i, total));}
            checkpoint("quadratic init");
        }
        // Read a product out of its transformed branches, applying the inverse-transform scale.
        // For d = 1 the branches hold the product modulo x^n - i and modulo x^n + i, so the low
        // half of the result is their half-sum and the high half their half-difference over i.
        // The half of the result that the branches are recombined into is parked in the output
        // itself, which is always long enough: the high half is only wanted where the low half
        // has already been read back.
        static void recover(std::array<cvector, 2> const& parts, size_t n, double factor, auto& out, size_t k) {
            // Residues are raw 32-bit words exactly when the modulus is a compile-time constant,
            // which is what lets the recombination stay vectorized.
            constexpr bool plain = fixed_mod && std::ranges::contiguous_range<std::decay_t<decltype(out)>>;
            constexpr bool raw = plain && sizeof(base) == 4;
            constexpr bool wide = plain && sizeof(base) == 8;
            const lattice L;
            auto scale = vz + factor;
            size_t low = std::min(k, n);
            auto project_at = [&](cvector const& part, size_t i, bool negative) {
                if(d == 1) {return project<true>(part.at(i) * scale, negative, L);}
                else {return project<false>(part.at(i) * scale, negative, L);}
            };
            // Eight coefficients at a time, the width the recombination works in.
            auto project8 = [&](cvector const& part, size_t i, bool negative, size_t count) {
                auto lo8 = project_at(part, i, negative);
                auto hi8 = count > flen ? project_at(part, i + flen, negative) : u32x4{};
                return __builtin_shufflevector(lo8, hi8, 0, 1, 2, 3, 4, 5, 6, 7);
            };
            auto store8 = [&](size_t idx, u32x8 v, size_t count) {
                if constexpr(raw) {
                    if(count == 8) {std::memcpy(reinterpret_cast<uint32_t*>(std::data(out)) + idx, &v, sizeof(v)); return;}
                } else if constexpr(wide) {
                    if(count == 8) {
                        auto words = reinterpret_cast<uint64_t*>(std::data(out)) + idx;
                        auto lo4 = __builtin_convertvector(u32x4{v[0], v[1], v[2], v[3]}, u64x4);
                        auto hi4 = __builtin_convertvector(u32x4{v[4], v[5], v[6], v[7]}, u64x4);
                        std::memcpy(words, &lo4, sizeof(lo4));
                        std::memcpy(words + flen, &hi4, sizeof(hi4));
                        return;
                    }
                }
                for(size_t l = 0; l < count; l++) {out[idx + l].setr(typename base::UInt(v[l]));}
            };
            auto load8 = [&](size_t idx, size_t count) {
                u32x8 v{};
                if constexpr(raw) {
                    if(count == 8) {std::memcpy(&v, reinterpret_cast<uint32_t const*>(std::data(out)) + idx, sizeof(v)); return v;}
                } else if constexpr(wide) {
                    if(count == 8) {
                        auto words = reinterpret_cast<uint64_t const*>(std::data(out)) + idx;
                        u64x4 lo4, hi4;
                        std::memcpy(&lo4, words, sizeof(lo4));
                        std::memcpy(&hi4, words + flen, sizeof(hi4));
                        auto lo = __builtin_convertvector(lo4, u32x4), hi = __builtin_convertvector(hi4, u32x4);
                        return u32x8(__builtin_shufflevector(lo, hi, 0, 1, 2, 3, 4, 5, 6, 7));
                    }
                }
                for(size_t l = 0; l < count; l++) {v[l] = uint32_t(out[idx + l].getr());}
                return v;
            };
            for(size_t i = 0; i < low; i += 8) {
                store8(i, project8(parts[0], i, false, std::min<size_t>(8, low - i)), std::min<size_t>(8, low - i));
            }
            if(d > 1) {checkpoint("quadratic recover"); return;}
            const uint32_t mod32 = uint32_t(base::mod()), imod32 = -inv2<uint32_t>(base::mod());
            auto highmul = u32x8{} + uint32_t(((base(2) * base(root)).inv() * bpow(base(2), 32)).getr());
            for(size_t i = 0; i < low; i += 8) {
                size_t count = std::min<size_t>(8, low - i);
                auto minus = project8(parts[1], i, true, count);
                auto plus = load8(i, count);
                auto lo8 = reduce_once(plus + minus, mod32);
                lo8 = (lo8 + (lo8 & 1) * mod32) >> 1;
                auto hi8 = reduce_once(montgomery_mul(plus + mod32 - minus, highmul, mod32, imod32), mod32);
                store8(i, lo8, count);
                if(n + i < k) {store8(n + i, hi8, std::min<size_t>(count, k - n - i));}
            }
            checkpoint("quadratic recover");
        }
        // Cyclic product modulo x^k - 1, in place over a, for operands of exactly k coefficients.
        // The transform evaluates at the k-th roots of i, so twisting coefficient j by w^j with
        // w^k = -i turns the product it computes, modulo x^k - i, into the cyclic one; the
        // inverse twist is folded into the readback. A fold of the linear product would need
        // twice the transform.
        static void cyclic(auto& a, auto const& b, size_t k) {
            init();
            assert(available && std::popcount(k) == 1 && std::size(a) == k && std::size(b) == k);
            bool square = (void const*)std::data(a) == (void const*)std::data(b);
            // w^j from a table of every fourth power, built the way the root tables of the
            // transform are, and four consecutive powers per vector by one broadcast multiply.
            size_t groups = k / flen;
            size_t fine_bits = std::min<size_t>(8, std::countr_zero(std::max<size_t>(groups, 1)));
            size_t fine = size_t(1) << fine_bits;
            big_vector<point> low(fine), high((groups + fine - 1) >> fine_bits), step(groups);
            auto w = [&](size_t j) {return polar<ftype>(1., -std::numbers::pi * ftype(j) / ftype(2 * k));};
            for(size_t t = 0; t < low.size(); t++) {low[t] = w(flen * t);}
            for(size_t c = 0; c < high.size(); c++) {high[c] = w(flen * (c << fine_bits));}
            for(size_t c = 0; c < groups; c++) {step[c] = high[c >> fine_bits] * low[c & (fine - 1)];}
            vpoint quarter, quarter_conj;
            for(size_t l = 0; l < flen; l++) {
                point v = w(l);
                real(quarter)[l] = real(v); imag(quarter)[l] = imag(v);
                real(quarter_conj)[l] = real(v); imag(quarter_conj)[l] = -imag(v);
            }
            auto twiddle = [&](size_t j, bool inverse, ftype scale) {
                point t = step[j / flen];
                vpoint head = {vz + real(t) * scale, vz + (inverse ? -imag(t) : imag(t)) * scale};
                return head * (inverse ? quarter_conj : quarter);
            };
            auto run = [&]<bool unit>() {
                const lattice L;
                cvector A(k), B(0);
                lift<unit>(A, a, k, false, seed());
                if(!square) {lift<unit>(B, b, k, false, seed());}
                for(size_t j = 0; j < k; j += flen) {
                    auto w = twiddle(j, false, 1);
                    A.at(j) *= w;
                    if(!square) {B.at(j) *= w;}
                }
                A.forward();
                if(!square) {B.forward();}
                A.multiply(square ? A : B);
                for(size_t j = 0; j < k; j += flen) {
                    auto v = project<unit>(A.at(j) * twiddle(j, true, ftype(flen) / ftype(k)), false, L);
                    for(size_t l = 0; l < flen; l++) {a[j + l].setr(typename base::UInt(v[l]));}
                }
            };
            if(d == 1) {run.template operator()<true>();} else {run.template operator()<false>();}
            checkpoint("quadratic recover");
        }
        static u64x4 seed() {
            return u64x4{random::rng() | 1, random::rng() | 1, random::rng() | 1, random::rng() | 1};
        }
        // a <- a * b for d = 1, or a <- a * a with square set (b is not read then).
        // The first branch is stored in the upper half of a, so the part of a that wraps around
        // is set aside first.
        static void mul_branches(auto& a, auto const& b, bool square) {
            size_t as = std::size(a), bs = square ? as : std::size(b), need = as + bs - 1;
            size_t n = length(as, bs);
            assert(available && d == 1);
            // Coefficients from 2n on come back negated at the bottom; there are few of them.
            std::array<uint32_t, short_tail> high{};
            size_t tail = need > 2 * n ? need - 2 * n : 0;
            for(size_t i = 0; i < tail; i++) {
                auto const* x = reinterpret_cast<const uint32_t*>(std::data(a));
                auto const* y = square ? x : reinterpret_cast<const uint32_t*>(std::data(b));
                uint64_t sum = 0;
                for(size_t j = 2 * n + i - bs + 1; j < as; j++) {
                    sum = (sum + uint64_t(x[j]) * y[2 * n + i - j]) % prime;
                }
                high[i] = uint32_t(sum);
            }
            // Montgomery constants for the recombination of the two branches.
            const uint32_t mod32 = uint32_t(base::mod()), imod32 = -inv2<uint32_t>(base::mod());
            base r32 = bpow(base(2), 32);
            auto highmul = u32x8{} + uint32_t(((base(2) * base(root)).inv() * r32).getr());
            u64x4 seed_a = seed(), seed_b = seed();
            const lattice L;
            big_vector<base> a_upper(std::begin(a) + std::min(as, n), std::end(a));
            a.resize(2 * n);
            std::span<base const> a_lower = std::span(a).first(std::min(as, n));
            // A square may come with b aliasing the storage that a has just left.
            auto b_lower = square ? a_lower : std::span<base const>(b).first(std::min(bs, n));
            auto b_upper = square ? std::span<base const>(a_upper) : std::span<base const>(b).subspan(b_lower.size());
            cvector A(0), B(0);
            auto* out = reinterpret_cast<uint32_t*>(std::data(a));
            for(bool negative: {false, true}) {
                if(n == (1 << 24)) {
                    if constexpr(cvector::fuse_forward) {
                        cvector::fuse_args fa{out, as, seed_a[0], L.a, L.b, L.a / double(base::mod()), L.b / double(base::mod())};
                        cvector::fuse_args fb{reinterpret_cast<const uint32_t*>(std::data(b)), bs, seed_b[0], fa.a, fa.b, fa.a_over_p, fa.b_over_p};
                        if(negative) {A.template cache_product<true>(B, fa, fb);}
                        else {A.template cache_product<false>(B, fa, fb);}
                    } else {
                        // cache_product transforms both operands in place, so a square is lifted twice.
                        fill<true, true>(A, a_lower, a_upper, n, negative, seed_a);
                        fill<true, true>(B, b_lower, b_upper, n, negative, seed_b);
                        if(negative) {A.template cache_product<true>(B);}
                        else {A.template cache_product<false>(B);}
                    }
                } else {
                    fill<false, true>(A, a_lower, a_upper, n, negative, seed_a);
                    if(!square) {fill<false, true>(B, b_lower, b_upper, n, negative, seed_b);}
                    A.multiply(square ? A : B);
                }
                auto scale = vz + double(flen) / double(n);
                for(size_t i = 0; i < n; i += 8) {
                    auto sum0 = project<true>(A.at(i) * scale, negative, L);
                    auto sum1 = project<true>(A.at(i + 4) * scale, negative, L);
                    auto sum = __builtin_shufflevector(sum0, sum1, 0, 1, 2, 3, 4, 5, 6, 7);
                    if(negative) {
                        u32x8 plus;
                        std::memcpy(&plus, out + n + i, sizeof(plus));
                        auto lo = plus + sum;
                        // Unsigned reduction: these sums pass 2^31 for moduli above 2^30.
                        lo = reduce_once(lo, mod32);
                        lo = (lo + (lo & 1) * base::mod()) >> 1;
                        auto hi = reduce_once(montgomery_mul(plus + base::mod() - sum, highmul, mod32, imod32), mod32);
                        std::memcpy(out + i, &lo, sizeof(lo));
                        std::memcpy(out + n + i, &hi, sizeof(hi));
                    } else {
                        std::memcpy(out + n + i, &sum, sizeof(sum));
                    }
                }
                checkpoint("quadratic recover");
            }
            a.resize(need);
            out = reinterpret_cast<uint32_t*>(std::data(a));
            for(size_t i = 0; i < tail; i++) {
                uint32_t sum = out[i] + high[i];
                out[i] = std::min(sum, sum - prime);
                out[2 * n + i] = high[i];
            }
        }
        // a <- a * b for d > 1, or a <- a * a with square set: i is not in the ring, so two
        // branches could only be recombined as complex numbers, which is slower and less exact
        // than one full transform.
        static void mul_single(auto& a, auto const& b, bool square) {
            size_t as = std::size(a), bs = square ? as : std::size(b), need = as + bs - 1;
            size_t n = length(as, bs), tail = need > n ? need - n : 0;
            assert(available && d > 1);
            const lattice L;
            cvector A(0), B(0);
            std::span<base const> none;
            // Coefficients from n on, in ring coordinates: they come back multiplied by i.
            std::array<point, short_tail> high{};
            auto wrapped = [&](cvector const& rhs) {
                for(size_t i = 0; i < tail; i++) {
                    for(size_t j = n + i - bs + 1; j < as; j++) {
                        high[i] += A.template get<point>(j) * rhs.template get<point>(n + i - j);
                    }
                }
            };
            if(n == (1 << 24)) {
                fill<true, false>(A, a, none, n, false, seed());
                if(square) {fill<true, false>(B, a, none, n, false, seed());}
                else {fill<true, false>(B, b, none, n, false, seed());}
                wrapped(B);
                A.template cache_product<false>(B);
            } else {
                fill<false, false>(A, a, none, n, false, seed(), !tail);
                if(!square) {fill<false, false>(B, b, none, n, false, seed(), !tail);}
                if(tail) {
                    wrapped(square ? A : B);
                    A.forward();
                    if(!square) {B.forward();}
                }
                A.multiply(square ? A : B);
            }
            for(size_t i = 0; i < tail; i++) {
                A.set(i, A.template get<point>(i) - point(0, 1) * high[i] * (double(n) / double(flen)));
            }
            a.resize(n);
            auto* out = reinterpret_cast<uint32_t*>(std::data(a));
            auto scale = vz + double(flen) / double(n);
            for(size_t i = 0; i < n; i += flen) {
                auto sum = project<false>(A.at(i) * scale, false, L);
                std::memcpy(out + i, &sum, sizeof(sum));
            }
            checkpoint("quadratic recover");
            a.resize(need);
            out = reinterpret_cast<uint32_t*>(std::data(a));
            for(size_t i = 0; i < tail; i += flen) {
                vpoint lanes = {vz, vz};
                for(size_t j = i; j < std::min(tail, i + flen); j++) {
                    real(lanes)[j - i] = real(high[j]);
                    imag(lanes)[j - i] = imag(high[j]);
                }
                auto sum = project<false>(lanes, false, L);
                for(size_t j = i; j < std::min(tail, i + flen); j++) {out[n + j] = sum[j - i];}
            }
        }
        // Both routines read and write the storage as plain residues, which it is for modint<m>.
        // A runtime-modulus type keeps another form, so its operands are converted on the way.
        static void mul(auto& a, auto const& b, bool square) {
            static_assert(sizeof(std::decay_t<decltype(a[0])>) == 4);
            auto run = [&](auto const& rhs) {
                if constexpr(fixed_d == 1) {mul_branches(a, rhs, square);}
                else if constexpr(fixed_d > 1) {mul_single(a, rhs, square);}
                else if(d == 1) {mul_branches(a, rhs, square);}
                else {mul_single(a, rhs, square);}
            };
            if constexpr(fixed_mod) {run(b);}
            else {
                big_vector<base> plain;
                if(!square) {plain.assign(std::begin(b), std::end(b));}
                for(auto& x: plain) {x.setr_direct(x.getr());}
                for(auto& x: a) {x.setr_direct(x.getr());}
                run(plain);
                for(auto& x: a) {x.setr(x.getr_direct());}
            }
        }
    };

    // A polynomial transformed for repeated multiplication in the ring of quadratic<base>.
    //
    // capacity() is the number of product coefficients the transform represents, and both
    // operands of a product must fit in it. For d = 1 the polynomial is held as two conjugate
    // branches of capacity/2 points, modulo x^n - i and modulo x^n + i, whose moduli multiply
    // to x^(2n) + 1; for larger d as one transform of capacity points modulo x^n - i. Either
    // way a product of at most capacity coefficients does not wrap around.
    //
    // complete selects the layout. The default stops one radix level short, which the fused
    // product of multiply() finishes; the complete transform makes a product pointwise, which
    // is what an accumulated sum of products needs.
    //
    // A transform costs two doubles per coefficient, so the type is move-only and every copy is
    // spelled clone(). multiply() and square() consume the object and only read the other
    // operand, which is what lets one transform serve several products.
    //
    // A transform belongs to the modulus it was built under; a runtime modulus must not change
    // while one is alive.
    template<modint_type base, bool complete = false>
    struct spectrum {
        using ring = quadratic<base>;
        spectrum(auto const& a, size_t capacity): cap(std::max(2 * flen, std::bit_ceil(capacity))) {
            ring::init();
            assert(ring::available && std::size(a) <= cap);
            auto state = ring::seed();
            for(size_t j = 0; j < branches(); j++) {
                if(ring::d == 1) {ring::template lift<true>(parts[j], a, points(), j == 1, state);}
                else {ring::template lift<false>(parts[j], a, points(), false, state);}
                if constexpr(complete) {parts[j].template fft<false>();}
                else {parts[j].forward();}
            }
        }
        spectrum(spectrum&&) = default;
        spectrum& operator=(spectrum&&) = default;
        spectrum clone() const {return *this;}
        size_t capacity() const {return cap;}

        // out[0..k) = *this * other. Consumes *this and only reads other.
        void multiply(spectrum const& other, auto& out, size_t k) && requires(!complete) {
            assert(other.cap == cap && k <= cap);
            for(size_t j = 0; j < branches(); j++) {parts[j].multiply(other.parts[j]);}
            ring::recover(parts, points(), double(flen) / double(points()), out, k);
        }
        // out[0..k) = *this * *this, consuming *this.
        void square(auto& out, size_t k) && requires(!complete) {
            std::move(*this).multiply(*this, out, k);
        }
    private:
        spectrum(spectrum const&) = default;
        size_t branches() const {return ring::d == 1 ? 2 : 1;}
        size_t points() const {return ring::d == 1 ? cap / 2 : cap;}
        template<modint_type> friend struct product;
        std::array<cvector, 2> parts = {cvector(0), cvector(0)};
        size_t cap;
    };

    // A sum of products of transforms, read back once. Every term costs one pointwise pass and
    // the sum one inverse transform, which is what makes it worth keeping the operands around.
    template<modint_type base>
    struct product {
        using ring = quadratic<base>;
        using operand = spectrum<base, true>;
        explicit product(size_t capacity): cap(std::max(2 * flen, std::bit_ceil(capacity))) {
            ring::init();
            assert(ring::available);
            for(size_t j = 0; j < branches(); j++) {acc[j] = cvector(points());}
        }
        // *this += x * y.
        void add(operand const& x, operand const& y) {
            assert(x.cap == cap && y.cap == cap);
            for(size_t j = 0; j < branches(); j++) {
                for(size_t i = 0; i < points(); i += flen) {
                    acc[j].at(i) += x.parts[j].at(i) * y.parts[j].at(i);
                }
            }
            checkpoint("dot");
        }
        // out[0..k) = the accumulated sum, consuming *this.
        void recover(auto& out, size_t k) && {
            assert(k <= cap);
            for(size_t j = 0; j < branches(); j++) {acc[j].template ifft<false>();}
            ring::recover(acc, points(), 1.0, out, k);
        }
    private:
        size_t branches() const {return ring::d == 1 ? 2 : 1;}
        size_t points() const {return ring::d == 1 ? cap / 2 : cap;}
        std::array<cvector, 2> acc = {cvector(0), cvector(0)};
        size_t cap;
    };
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_RING_HPP
