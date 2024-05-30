#ifndef CP_ALGO_MATH_FFT_HPP
#define CP_ALGO_MATH_FFT_HPP
#include "common.hpp"
#include "modint.hpp"
#include <algorithm>
#include <complex>
#include <cassert>
#include <vector>
#include <bit>
namespace cp_algo::math::fft {
    using ftype = double;
    using point = std::complex<ftype>;

    std::vector<point> w; // w[2^n + k] = exp(pi * k / (2^n))
    std::vector<int> bitr;// b[2^n + k] = bitreverse(k)
    const ftype pi = acos(-1);
    bool initiated = 0;
    size_t bitreverse(size_t n, size_t k) {
        size_t hn = n / 2;
        if(k >= hn && n > 1) {
            return 2 * bitr[k] + 1;
        } else {
            return 2 * bitr[hn + k];
        }
    }
    void init() {
        if(!initiated) {
            w.resize(maxn);
            bitr.resize(maxn);
            for(int i = 1; i < maxn; i *= 2) {
                int ti = i / 2;
                ftype arg = pi / i;
                point base = std::polar(1., arg);
                point cur = 1.;
                for(int j = 0; j < i; j++) {
                    if((j & 15) == 0) {
                        cur = std::polar(1., j * arg);
                    }
                    w[i + j] = cur;
                    cur *= base;
                    if(ti) {
                        bitr[i + j] = 2 * bitr[ti + j % ti] + (j >= ti);
                    }
                }
            }
            initiated = 1;
        }
    }
    
    void ifft(auto &a, int n) {
        init();
        if(n == 1) {
            return;
        }
        for(int i = 1; i < n; i *= 2) {
            for(int j = 0; j < n; j += 2 * i) {
                for(int k = j; k < j + i; k++) {
                    std::tie(a[k], a[k + i]) = std::pair{
                        a[k] + a[k + i] * conj(w[i + k - j]),
                        a[k] - a[k + i] * conj(w[i + k - j])
                    };
                }
            }
        }
    }
    void fft(auto &a, int n) {
        init();
        if(n == 1) {
            return;
        }
        for(int i = n / 2; i >= 1; i /= 2) {
            for(int j = 0; j < n; j += 2 * i) {
                for(int k = j; k < j + i; k++) {
                    std::tie(a[k], a[k + i]) = std::pair{
                        a[k] + a[k + i],
                        (a[k] - a[k + i]) * w[i + k - j]
                    };
                }
            }
        }
    }
    template<typename base>
    void mul_slow(std::vector<base> &a, const std::vector<base> &b) {
        if(a.empty() || b.empty()) {
            a.clear();
        } else {
            int n = a.size();
            int m = b.size();
            a.resize(n + m - 1);
            for(int k = n + m - 2; k >= 0; k--) {
                a[k] *= b[0];
                for(int j = std::max(k - n + 1, 1); j < std::min(k + 1, m); j++) {
                    a[k] += a[k - j] * b[j];
                }
            }
        }
    }
    
    template<typename base>
    struct dft {
        std::vector<point> A;
        
        dft(std::vector<base> const& a, size_t n): A(n) {
            for(size_t i = 0; i < std::min(n, a.size()); i++) {
                A[i] = a[i];
            }
            if(n) {
                fft(A, n);
            }
        }
    
        std::vector<base> operator *= (dft const& B) {
            assert(A.size() == B.A.size());
            size_t n = A.size();
            if(!n) {
                return std::vector<base>();
            }
            for(size_t i = 0; i < n; i++) {
                A[i] *= B[i];
            }
            ifft(A, n);
            for(size_t i = 0; i < n; i++) {
                A[i] /= n;
            }
            if constexpr (std::is_same_v<base, point>) {
                return A;
            } else {
                return {begin(A), end(A)};
            }
        }

        auto operator * (dft const& B) const {
            return dft(*this) *= B;
        }
        
        point& operator [](int i) {return A[i];}
        point operator [](int i) const {return A[i];}
    };

    template<modint_type base>
    struct dft<base> {
        static constexpr int split = 1 << 15;
        std::vector<point> A;
        
        dft(std::vector<base> const& a, size_t n): A(n) {
            for(size_t i = 0; i < std::min(n, a.size()); i++) {
                A[i] = point(a[i].rem() % split, a[i].rem() / split);
            }
            if(n) {
                fft(A, n);
            }
        }

        std::vector<base> operator *= (dft const& B) {
            assert(A.size() == B.A.size());
            size_t n = A.size();
            if(!n) {
                return std::vector<base>();
            }
            std::vector<point> C(n);
            for(size_t i = 0; 2 * i <= n; i++) {
                size_t j = (n - i) % n;
                size_t x = bitreverse(n, i);
                size_t y = bitreverse(n, j);
                std::tie(C[x], A[x], C[y], A[y]) = std::make_tuple(
                    A[x] * (B[x] + conj(B[y])),
                    A[x] * (B[x] - conj(B[y])),
                    A[y] * (B[y] + conj(B[x])),
                    A[y] * (B[y] - conj(B[x]))
                );
            }
            ifft(C, n);
            ifft(A, n);
            int t = 2 * n;
            std::vector<base> res(n);
            for(size_t i = 0; i < n; i++) {
                base A0 = llround(C[i].real() / t);
                base A1 = llround(C[i].imag() / t + A[i].imag() / t);
                base A2 = llround(A[i].real() / t);
                res[i] = A0 + A1 * split - A2 * split * split;
            }
            return res;
        }

        auto operator * (dft const& B) const {
            return dft(*this) *= B;
        }
        
        point& operator [](int i) {return A[i];}
        point operator [](int i) const {return A[i];}
    };
    
    size_t com_size(size_t as, size_t bs) {
        if(!as || !bs) {
            return 0;
        }
        return std::bit_ceil(as + bs - 1);
    }
    
    template<typename base>
    void mul(std::vector<base> &a, std::vector<base> const& b) {
        if(std::min(a.size(), b.size()) < magic) {
            mul_slow(a, b);
            return;
        }
        auto n = com_size(a.size(), b.size());
        auto A = dft<base>(a, n);
        if(a == b) {
            a = A *= A;
        } else {
            a = A *= dft<base>(b, n);
        }
    }
    template<typename base>
    void circular_mul(std::vector<base> &a, std::vector<base> const& b) {
        auto n = std::bit_ceil(a.size());
        auto A = dft<base>(a, n);
        if(a == b) {
            a = A *= A;
        } else {
            a = A *= dft<base>(b, n);
        }
    }
}
#endif // CP_ALGO_MATH_FFT_HPP
