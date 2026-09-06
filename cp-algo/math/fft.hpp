#ifndef CP_ALGO_MATH_FFT_HPP
#define CP_ALGO_MATH_FFT_HPP
#include "dft.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::fft {
    void mul_slow(auto &a, auto const& b, size_t k) {
        if(!std::empty(a) && std::data(a) == std::data(b)) {
            auto copy = big_vector<std::decay_t<decltype(b[0])>>(begin(b), end(b));
            return mul_slow(a, copy, k);
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
            auto B = dft<base>(b | std::views::take(k), n / 2);
            a.resize((k + flen - 1) / flen * flen);
            A.mul_inplace(B, a, std::min(k, n));
            auto wrap = bpow(dft<base>::factor, n);
            for(size_t i = 0; i < tail; i++) {
                a[i] += wrap * high[i];
                if(n + i < k) {a[n + i] = high[i];}
            }
            a.resize(k);
            return;
        }
        auto A = dft<base>(a | std::views::take(k), n);
        if(std::data(a) == std::data(b)) {
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
        if(k <= (1 << 16)) {
            big_vector<base> ap(begin(a), end(a));
            mul_truncate(ap, b, 2 * k);
            mod_split(ap, k, bpow(dft<base>::factor, k));
            std::ranges::copy(ap | std::views::take(k), begin(a));
            return;
        }
        k /= 2;
        auto factor = bpow(dft<base>::factor, k);
        if(zero_upper) {
            std::ranges::copy(std::span(a).first(k), begin(a) + k);
            std::ranges::copy(std::span(b).first(k), begin(b) + k);
        } else {
            mod_split(a, k, factor);
            mod_split(b, k, factor);
        }
        auto la = std::span(a).first(k);
        auto lb = std::span(b).first(k);
        auto ra = std::span(a).last(k);
        auto rb = std::span(b).last(k);
        cyclic_mul(la, lb, k);
        auto A = dft<base>(ra, k / 2);
        auto B = dft<base>(rb, k / 2);
        A.mul_inplace(B, ra, k);
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
    void mul(auto &a, auto &&b) {
        if(std::empty(a) || std::empty(b)) {a.clear(); return;}
        if(std::data(a) == std::data(b)) {
            auto copy = make_copy(b);
            return mul(a, copy);
        }
        size_t N = size(a) + size(b);
        if(N > (1 << 20)) {
            N--;
            size_t NN = std::bit_ceil(N);
            bool zero_upper = std::max(size(a), size(b)) <= NN / 2;
            a.resize(NN);
            b.resize(NN);
            cyclic_mul(a, b, NN, zero_upper);
            a.resize(N);
        } else {
            mul_truncate(a, b, N - 1);
        }
    }
    void mul(auto &a, auto const& b) {
        if(std::empty(a) || std::empty(b)) {a.clear(); return;}
        size_t N = size(a) + size(b);
        if(N > (1 << 20)) {
            mul(a, make_copy(b));
        } else {
            mul_truncate(a, b, N - 1);
        }
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_FFT_HPP
