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
        auto A = dft<base>(a | std::views::take(k), n);
        auto B = dft<base>(b | std::views::take(k), n);
        a.resize((k + flen - 1) / flen * flen);
        A.mul_inplace(B, a, k);
        a.resize(k);
    }

    // store mod x^n-k in first half, x^n+k in second half
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
            xr = montgomery_mul(xr, cur, dft<base>::mod, dft<base>::imod);
            xr = xr >= base::mod() ? xr - base::mod() : xr;
            auto t = xr;
            xr = xl - t;
            xl += t;
            xl = xl >= base::mod() ? xl - base::mod() : xl;
            xr = xr >= base::mod() ? xr + base::mod() : xr;
            for(size_t k = 0; k < flen; k++) {
                x[i + k].setr(typename base::UInt(xl[k]));
                x[n + i + k].setr(typename base::UInt(xr[k]));
            }
        }
        cp_algo::checkpoint("mod split");
    }
    void cyclic_mul(auto &a, auto &&b, size_t k) {
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
        mod_split(a, k, factor);
        mod_split(b, k, factor);
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
        for(size_t i = 0; i < k; i++) {
            auto t = (a[i] + a[i + k]) * i2;
            a[i + k] = (a[i] - a[i + k]) * factor;
            a[i] = t;
        }
        cp_algo::checkpoint("mod join");
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
            a.resize(NN);
            b.resize(NN);
            cyclic_mul(a, b, NN);
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
