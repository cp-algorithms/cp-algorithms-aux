#ifndef CP_ALGO_MATH_FFT_HPP
#define CP_ALGO_MATH_FFT_HPP
#include "dft.hpp"
#include <cstring>
#include <tuple>
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::fft {
    void mul_slow(auto &a, auto const& b, size_t k) {
        if(!std::empty(a) && std::data(a) == std::data(b)) {
            using base = std::decay_t<decltype(a[0])>;
            size_t n = std::min(k, std::size(a)), m = std::min(k, std::size(b));
            if(!m) {a.clear(); return;}
            a.resize(k);
            // Descending output only reads original coefficients at indices <=j.
            for(size_t j = k; j-- > 0;) {
                base sum = 0;
                size_t lo = j >= n ? j + 1 - n : 0, hi = std::min(j + 1, m);
                for(size_t i = lo; i < hi; i++) {
                    if(n == m && i > j - i) {break;}
                    auto term = a[i] * a[j - i];
                    sum += n == m && i != j - i ? term + term : term;
                }
                a[j] = sum;
            }
            return;
        }
        if(std::empty(a) || std::empty(b)) {
            a.clear();
        } else {
            size_t n = std::min(k, std::size(a));
            size_t m = std::min(k, std::size(b));
            a.resize(k);
            for(int j = int(k - 1); j >= 0; j--) {
                a[j] *= b[0];
                for(int i = std::max(j - (int)n, 0) + 1; i < std::min(j + 1, (int)m); i++) {
                    a[j] += a[j - i] * b[i];
                }
            }
        }
    }
    size_t com_size(size_t as, size_t bs) {
        if(!as || !bs) {
            return 0;
        }
        return std::max(flen, std::bit_ceil(as + bs - 1) / 2);
    }
    void mul_truncate(auto &a, auto const& b, size_t k) {
        using base = std::decay_t<decltype(a[0])>;
        if(std::min({k, std::size(a), std::size(b)}) < magic) {
            mul_slow(a, b, k);
            return;
        }
        auto n = std::max(flen, std::bit_ceil(
            std::min(k, std::size(a)) + std::min(k, std::size(b)) - 1
        ) / 2);
        size_t as = std::min(k, std::size(a)), bs = std::min(k, std::size(b));
        size_t tail = as + bs - 1 - n;
        // Correct a short wrapped tail instead of doubling the FFT size.
        if(tail <= 32 && as <= n && bs <= n) {
            std::array<base, 32> high{};
            for(size_t i = 0; i < tail; i++) {
                for(size_t j = n + i - bs + 1; j < as; j++) {
                    high[i] += a[j] * b[n + i - j];
                }
            }
            auto A = dft<base>(a | std::views::take(k), n / 2);
            if(as == bs && std::data(a) == std::data(b)) {
                a.resize((k + flen - 1) / flen * flen);
                A.mul(A, a, std::min(k, n));
            } else {
                auto B = dft<base>(b | std::views::take(k), n / 2);
                a.resize((k + flen - 1) / flen * flen);
                A.mul_inplace(B, a, std::min(k, n));
            }
            auto wrap = bpow(dft<base>::factor, n);
            for(size_t i = 0; i < tail; i++) {
                a[i] += wrap * high[i];
                if(n + i < k) {a[n + i] = high[i];}
            }
            a.resize(k);
            return;
        }
        auto A = dft<base>(a | std::views::take(k), n);
        if(as == bs && std::data(a) == std::data(b)) {
            a.resize((k + flen - 1) / flen * flen);
            A.mul(A, a, k);
        } else {
            auto B = dft<base>(b | std::views::take(k), n);
            a.resize((k + flen - 1) / flen * flen);
            A.mul_inplace(B, a, k);
        }
        a.resize(k);
    }

    // store mod x^n-k in first half, x^n+k in second half
    // inverse reconstructs the halves with k = 1/(2 * forward_k).
    template<bool inverse = false>
    void mod_split(auto &&x, size_t n, auto k) {
        using base = std::decay_t<decltype(k)>;
        dft<base>::init();
        assert(std::size(x) == 2 * n);
        u64x4 cur = u64x4{} + (k * bpow(base(2), 32)).getr();
        for(size_t i = 0; i < n; i += flen) {
            u64x4 xl = {
                x[i].getr(),
                x[i + 1].getr(),
                x[i + 2].getr(),
                x[i + 3].getr()
            };
            u64x4 xr = {
                x[n + i].getr(),
                x[n + i + 1].getr(),
                x[n + i + 2].getr(),
                x[n + i + 3].getr()
            };
            if constexpr(!inverse) {
                xr = montgomery_mul(xr, cur, dft<base>::mod, dft<base>::imod);
                xr = xr >= base::mod() ? xr - base::mod() : xr;
            }
            auto t = xr;
            xr = xl - t;
            xl += t;
            xl = xl >= base::mod() ? xl - base::mod() : xl;
            xr = xr >= base::mod() ? xr + base::mod() : xr;
            if constexpr(inverse) {
                xl = (xl + (xl & 1) * base::mod()) >> 1;
                xr = montgomery_mul(xr, cur, dft<base>::mod, dft<base>::imod);
                xr = xr >= base::mod() ? xr - base::mod() : xr;
            }
            for(size_t k = 0; k < flen; k++) {
                x[i + k].setr(typename base::UInt(xl[k]));
                x[n + i + k].setr(typename base::UInt(xr[k]));
            }
        }
        cp_algo::checkpoint(inverse ? "mod join" : "mod split");
    }
    // zero_upper skips arithmetic on the known zero padding in the first split.
    void cyclic_mul(auto &a, auto &&b, size_t k, bool zero_upper = false) {
        assert(std::popcount(k) == 1);
        assert(std::size(a) == std::size(b) && std::size(a) == k);
        using base = std::decay_t<decltype(a[0])>;
        dft<base>::init();
        bool square = std::data(a) == std::data(b);
        if(k <= (1 << 16)) {
            big_vector<base> ap(begin(a), end(a));
            if(square) {mul_truncate(ap, ap, 2 * k);}
            else {mul_truncate(ap, b, 2 * k);}
            mod_split(ap, k, bpow(dft<base>::factor, k));
            std::ranges::copy(ap | std::views::take(k), begin(a));
            return;
        }
        k /= 2;
        auto factor = bpow(dft<base>::factor, k);
        if(zero_upper) {
            std::ranges::copy(std::span(a).first(k), begin(a) + k);
            if(!square) {std::ranges::copy(std::span(b).first(k), begin(b) + k);}
        } else {
            mod_split(a, k, factor);
            if(!square) {mod_split(b, k, factor);}
        }
        auto la = std::span(a).first(k);
        auto lb = std::span(b).first(k);
        auto ra = std::span(a).last(k);
        auto rb = std::span(b).last(k);
        cyclic_mul(la, lb, k);
        auto A = dft<base>(ra, k / 2);
        if(square) {A.mul(A, ra, k);}
        else {
            auto B = dft<base>(rb, k / 2);
            A.mul_inplace(B, ra, k);
        }
        base i2 = base(2).inv();
        factor = factor.inv() * i2;
        mod_split<true>(a, k, factor);
    }
    auto make_copy(auto &&x) {
        return x;
    }
    void cyclic_mul(auto &a, auto const& b, size_t k) {
        return cyclic_mul(a, make_copy(b), k);
    }
    // Gaussian-integer convolution for a 32-bit prime p = a^2 + b^2 (p = 1 mod 4).
    // Residues become re + im*i with |re|, |im| < sqrt(p) after reducing by the lattice
    // spanned by (a, b) and (-b, a); i is identified with root = b / a, a square root of -1.
    // The product is one complex convolution instead of three real ones, computed modulo
    // x^n - i and x^n + i (the latter through conjugation) and recombined.
    template<modint_type base>
    struct gaussian {
        static inline bool ready = false, available = false;
        static inline uint32_t a = 0, b = 0, root = 0;
        static void init() {
            if(ready) {return;}
            ready = true;
            uint64_t p = base::mod();
            if(p % 4 != 1 || p >= (uint64_t(1) << 31)) {return;}
            // A square root of -1 from a quadratic non-residue, then Euclid down to sqrt(p).
            uint64_t s = 0;
            for(uint64_t g = 2; g < p; g++) {
                if(bpow(base(g), (p - 1) / 2) != base(1)) {s = bpow(base(g), (p - 1) / 4).getr(); break;}
            }
            uint64_t r0 = p, r1 = s;
            while(r1 * r1 > p) {std::tie(r0, r1) = std::pair{r1, r0 % r1};}
            uint64_t aa = r1, bb = uint64_t(std::sqrt((long double)(p - aa * aa)) + 0.5);
            if(aa * aa + bb * bb != p) {return;}
            a = uint32_t(aa); b = uint32_t(bb);
            root = (base(b) / base(a)).getr();
            assert(base(root) * base(root) == base(p - 1));
            available = true;
        }
        // Whether the Gaussian path is used for operands of these sizes.
        static bool usable(size_t as, size_t bs) {
            init();
            if(!available || !as || !bs) {return false;}
            size_t n = std::max(flen, std::bit_ceil(as + bs - 1) / 2);
            return std::max(as, bs) <= n && n >= (size_t(1) << 18);
        }
        // Map a rounded Gaussian coordinate pair back to the residue re + root*im.
        static u32x4 project(vpoint value, bool negative) {
            const double p = base::mod(), da = a, db = b, dr = root;
            auto R = round(real(value)), I = round(negative ? -imag(value) : imag(value));
            auto q = round(I * (1.0 / db));
            auto U = vftype(_mm256_fnmadd_pd(__m256d(q), _mm256_set1_pd(da), __m256d(R)));
            auto V = vftype(_mm256_fnmadd_pd(__m256d(q), _mm256_set1_pd(db), __m256d(I)));
            auto h = vftype(_mm256_fmadd_pd(__m256d(V), _mm256_set1_pd(dr), __m256d(U)));
            q = round(h * (1.0 / p));
            auto out = vftype(_mm256_fnmadd_pd(__m256d(q), _mm256_set1_pd(p), __m256d(h)));
            out = out < 0 ? out + p : out;
            return u32x4(_mm256_cvttpd_epi32(__m256d(out)));
        }
        // Lift residues to Gaussian coordinates with stochastic rounding (generic sizes).
        static void fill(cvector& c, auto const& x, size_t n, bool negative, uint64_t seed) {
            cvector::fuse_args fa{reinterpret_cast<const uint32_t*>(std::data(x)), std::size(x), seed,
                                  double(a), double(b), a / double(base::mod()), b / double(base::mod())};
            c.r.clear(); c.r.reserve(n / flen);
            for(size_t i = 0; i < std::size(x); i += flen) {
                i32x4 bits{};
                if(i + flen <= std::size(x)) {std::memcpy(&bits, fa.src + i, sizeof(bits));}
                else {for(size_t j = i; j < std::size(x); j++) {bits[j - i] = int32_t(fa.src[j]);}}
                auto v = __builtin_convertvector(bits, vftype);
                u32x4 h = u32x4{uint32_t(i), uint32_t(i + 1), uint32_t(i + 2), uint32_t(i + 3)} ^ uint32_t(seed);
                h *= 0x9E3779B1u; h ^= h >> 15; h *= 0x85EBCA77u; h ^= h >> 13; h *= 0xC2B2AE3Du; h ^= h >> 16;
                auto noise = __builtin_convertvector(i32x4(h), vftype) * 0x1p-32;
                auto q = round(v * fa.a_over_p + noise), t = round(v * fa.b_over_p + noise);
                auto re = v - q * fa.a - t * fa.b, im = t * fa.a - q * fa.b;
                c.r.push_back(vpoint{re, negative ? -im : im});
            }
            size_t old = c.r.size();
            c.r.resize(n / flen);
            std::fill(c.r.begin() + old, c.r.end(), vpoint{});
            checkpoint("gaussian init");
            if(n != (1 << 24)) {c.fft();}
        }
        // a <- a * b; both operands must fit in half of the padded product length.
        static void mul(auto& a, auto const& b) {
            static_assert(sizeof(std::decay_t<decltype(a[0])>) == 4);
            size_t as = std::size(a), bs = std::size(b), need = as + bs - 1;
            size_t n = std::max(flen, std::bit_ceil(need) / 2);
            assert(available && as <= n && bs <= n);
            using D = dft<base>;
            D::init();
            base r32 = bpow(base(2), 32);
            auto highmul = u32x8{} + uint32_t(((base(2) * base(root)).inv() * r32).getr());
            uint64_t seed_a = random::rng() | 1, seed_b = random::rng() | 1;
            a.resize(2 * n);
            cvector A(0), B(0);
            auto* out = reinterpret_cast<uint32_t*>(std::data(a));
            for(bool negative: {false, true}) {
                if(n == (1 << 24)) {
                    if constexpr(cvector::fuse_forward) {
                        cvector::fuse_args fa{out, as, seed_a, double(gaussian::a), double(gaussian::b),
                                              gaussian::a / double(base::mod()), gaussian::b / double(base::mod())};
                        cvector::fuse_args fb{reinterpret_cast<const uint32_t*>(std::data(b)), bs, seed_b, fa.a, fa.b, fa.a_over_p, fa.b_over_p};
                        if(negative) {A.template cache_product<true>(B, fa, fb);}
                        else {A.template cache_product<false>(B, fa, fb);}
                    } else {
                        fill(A, std::span(a).first(as), n, negative, seed_a);
                        fill(B, b, n, negative, seed_b);
                        if(negative) {A.template cache_product<true>(B);}
                        else {A.template cache_product<false>(B);}
                    }
                } else {
                    fill(A, std::span(a).first(as), n, negative, seed_a);
                    fill(B, b, n, negative, seed_b);
                    A.dot(B);
                    A.template ifft<true, false>();
                }
                using i32x8 = simd<int32_t, 8>;
                auto scale = vz + double(flen) / double(n);
                for(size_t i = 0; i < n; i += 8) {
                    auto sum0 = project(A.at(i) * scale, negative);
                    auto sum1 = project(A.at(i + 4) * scale, negative);
                    auto sum = __builtin_shufflevector(sum0, sum1, 0, 1, 2, 3, 4, 5, 6, 7);
                    if(negative) {
                        u32x8 plus;
                        std::memcpy(&plus, out + n + i, sizeof(plus));
                        auto lo = plus + sum;
                        lo = i32x8(lo) >= int32_t(base::mod()) ? lo - base::mod() : lo;
                        lo = (lo + (lo & 1) * base::mod()) >> 1;
                        auto hi = montgomery_mul(plus + base::mod() - sum, highmul, D::mod, D::imod);
                        hi = i32x8(hi) >= int32_t(base::mod()) ? hi - base::mod() : hi;
                        std::memcpy(out + i, &lo, sizeof(lo));
                        std::memcpy(out + n + i, &hi, sizeof(hi));
                    } else {
                        std::memcpy(out + n + i, &sum, sizeof(sum));
                    }
                }
                checkpoint("gaussian recover");
            }
            a.resize(need);
        }
    };
    namespace impl {
        // Overlap-add for a short fixed operand; every block reuses its transform.
        void mul_unbalanced(auto &a, auto const& b) {
            using base = std::decay_t<decltype(a[0])>;
            auto x = std::span<base const>(a), y = std::span<base const>(b);
            if(x.size() < y.size()) {std::swap(x, y);}
            constexpr size_t length = 1 << 15;
            size_t step = length - y.size() + 1;
            auto fixed = dft<base>(y, length / 2);
            std::decay_t<decltype(a)> result(x.size() + y.size() - 1);
            big_vector<base> work(length);
            for(size_t start = 0; start < x.size(); start += step) {
                size_t count = std::min(step, x.size() - start);
                auto block = dft<base>(x.subspan(start, count), length / 2);
                size_t need = count + y.size() - 1;
                block.mul(fixed, work, need);
                for(size_t i = 0; i < need; i++) {result[start + i] += work[i];}
            }
            a = std::move(result);
        }
    }
    void mul(auto &a, auto &&b) {
        if(std::empty(a) || std::empty(b)) {a.clear(); return;}
        bool square = std::data(a) == std::data(b) && std::size(a) == std::size(b);
        if(!square && std::data(a) == std::data(b)) {
            auto copy = make_copy(b);
            return mul(a, copy);
        }
        using base = std::decay_t<decltype(a[0])>;
        if constexpr(sizeof(base) == 4) {
            if(gaussian<base>::usable(std::size(a), std::size(b))) {
                if(square) {auto copy = make_copy(b); return gaussian<base>::mul(a, copy);}
                return gaussian<base>::mul(a, b);
            }
        }
        size_t small = std::min(size(a), size(b)), large = std::max(size(a), size(b));
        if(small >= magic && small <= 4096 && large >= (1 << 20) && large / small >= 64) {
            return impl::mul_unbalanced(a, b);
        }
        size_t N = size(a) + size(b);
        if(N > (1 << 20)) {
            N--;
            size_t NN = std::bit_ceil(N);
            bool zero_upper = std::max(size(a), size(b)) <= NN / 2;
            a.resize(NN);
            // Compute the negative branch before the positive branch consumes the inputs.
            // Only the result needs the upper half; b never needs duplicated padding.
            if(zero_upper && !square) {
                size_t half = NN / 2;
                b.resize(half);
                auto lo = std::span(a).first(half), hi = std::span(a).last(half);
                {
                    auto A = dft<base>(lo, half / 2);
                    auto B = dft<base>(b, half / 2);
                    A.mul_inplace(B, hi, half);
                }
                cyclic_mul(lo, b, half);
                mod_split<true>(a, half, (base(2) * bpow(dft<base>::factor, half)).inv());
            } else {
                if(!square) {b.resize(NN);}
                cyclic_mul(a, b, NN, zero_upper);
            }
            a.resize(N);
        } else {
            mul_truncate(a, b, N - 1);
        }
    }
    void mul(auto &a, auto const& b) {
        if(std::empty(a) || std::empty(b)) {a.clear(); return;}
        size_t small = std::min(size(a), size(b)), large = std::max(size(a), size(b));
        if(small >= magic && small <= 4096 && large >= (1 << 20) && large / small >= 64) {
            return impl::mul_unbalanced(a, b);
        }
        using base = std::decay_t<decltype(a[0])>;
        if constexpr(sizeof(base) == 4) {
            if(gaussian<base>::usable(std::size(a), std::size(b))) {
                if(std::data(a) == std::data(b)) {auto copy = make_copy(b); return gaussian<base>::mul(a, copy);}
                return gaussian<base>::mul(a, b);
            }
        }
        size_t N = size(a) + size(b);
        if(N > (1 << 20)) {
            if(std::data(a) == std::data(b) && std::size(a) == std::size(b)) {mul(a, a);}
            else {mul(a, make_copy(b));}
        } else {
            mul_truncate(a, b, N - 1);
        }
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_FFT_HPP
