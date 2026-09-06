#include "cp-algo/math/fft.hpp"
#include <random>
#include <iostream>
using namespace cp_algo;
using namespace cp_algo::math;

template<int mod> void convolution_boundaries() {
    using T = modint<mod>;
    std::mt19937 rng(42);
    // In-place small squares and shorter aliased prefixes, including truncation.
    for(size_t n = 1; n <= 32; n++) {
        big_vector<T> a(n);
        for(auto &v: a) {v = rng() % mod;}
        for(size_t m = 0; m <= n; m++) {
            big_vector<T> c(n + m);
            for(size_t i = 0; i < n; i++) {
                for(size_t j = 0; j < m; j++) {c[i + j] += a[i] * a[j];}
            }
            for(size_t k = 0; k <= n + m + 2; k++) {
                auto x = a, want = c;
                want.resize(m ? k : 0);
                fft::mul_truncate(x, std::span(x).first(m), k);
                assert(x == want);
                if(m == n) {
                    x = a;
                    fft::mul_truncate(x, x, k);
                    assert(x == want);
                }
            }
        }
    }
    for(int n: {64, 65, 66, 79, 80, 81, 97, 127, 128, 129, 130, 145, 257, 513}) {
        for(int m: {64, 65, 66, 79, 97, 127, 128, 129, 257}) {
            big_vector<T> a(n), b(m), c(n + m - 1);
            for(auto &v: a) {v = rng() % mod;}
            for(auto &v: b) {v = rng() % mod;}
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < m; j++) {c[i + j] += a[i] * b[j];}
            }
            for(size_t k: {size_t(0), size_t(63), size_t(64), size_t(n), c.size() - 1, c.size(), c.size() + 7}) {
                auto x = a, want = c;
                want.resize(k);
                fft::mul_truncate(x, b, k);
                assert(x == want);
            }
        }
        big_vector<T> a(n), c(2 * n - 1);
        for(auto &v: a) {v = rng() % mod;}
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {c[i + j] += a[i] * a[j];}
        }
        for(size_t k: {size_t(n), c.size(), c.size() + 7}) {
            auto x = a, want = c;
            want.resize(k);
            fft::mul_truncate(x, x, k);
            assert(x == want);
        }
        // Sharing the first coefficient does not make a shorter prefix a square.
        for(size_t m: {size_t(31), size_t(64), size_t(n - 1)}) {
            m = std::min(m, a.size());
            big_vector<T> product(a.size() + m - 1);
            for(size_t i = 0; i < a.size(); i++) {
                for(size_t j = 0; j < m; j++) {product[i + j] += a[i] * a[j];}
            }
            for(size_t k: {size_t(63), size_t(n), product.size(), product.size() + 7}) {
                auto x = a, want = product;
                want.resize(k);
                fft::mul_truncate(x, std::span(x).first(m), k);
                assert(x == want);
            }
        }
    }
}

int main() {
    convolution_boundaries<998244353>();
    convolution_boundaries<1000000007>();
    convolution_boundaries<17>();
    convolution_boundaries<65537>();
    std::mt19937 rng(42);
    // Exercise both FFT layouts as the root cache grows, then reuse smaller sizes.
    for(size_t n: {4, 8, 16, 1 << 16, 1 << 20, 1 << 22, 8}) {
        auto roundtrip = [&]<bool partial, bool normalize>() {
            fft::cvector a(n);
            for(size_t i = 0; i < n; i++) {
                a.set(i, fft::point(int(rng() % 1001) - 500, int(rng() % 1001) - 500));
            }
            auto expected = a;
            a.fft<partial>();
            a.ifft<partial, normalize>();
            double scale = normalize ? 1 : double(partial ? n / fft::flen : n);
            for(size_t i = 0; i < n; i++) {
                assert(abs(a.get(i) / scale - expected.get(i)) < 1e-7);
            }
        };
        roundtrip.template operator()<true, true>();
        roundtrip.template operator()<false, true>();
        roundtrip.template operator()<true, false>();
        roundtrip.template operator()<false, false>();
    }
    auto split_roundtrip = [&]<int mod>() {
        using T = modint<mod>;
        big_vector<T> a(4096);
        for(size_t i = 0; i < a.size(); i++) {
            a[i] = i % 3 == 0 ? 0 : i % 3 == 1 ? mod - 1 : rng() % mod;
        }
        auto expected = a;
        fft::mod_split(a, a.size() / 2, T(17));
        fft::mod_split<true>(a, a.size() / 2, T(34).inv());
        assert(a == expected);
    };
    split_roundtrip.template operator()<998244353>();
    split_roundtrip.template operator()<1000000007>();

    // Multiplication by (1-x)^2 gives an exact linear-time reference above the large cutoff.
    using T = modint<998244353>;
    for(size_t n: {(1 << 20) - 1, (1 << 20) + 7}) {
        big_vector<T> a(n), b{1, -2, 1};
        for(auto &x: a) {x = rng() % T::mod();}
        auto original = a;
        fft::mul(a, b);
        for(size_t i = 0; i < a.size(); i++) {
            T expected = i < original.size() ? original[i] : T(0);
            if(i && i - 1 < original.size()) {expected -= original[i - 1] * 2;}
            if(i > 1 && i - 2 < original.size()) {expected += original[i - 2];}
            assert(a[i] == expected);
        }
    }
    std::cout << "FFT and large convolution properties passed\n";
}
