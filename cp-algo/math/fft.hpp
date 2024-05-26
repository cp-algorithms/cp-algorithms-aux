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
    void init() {
        if(!initiated) {
            w.resize(maxn);
            bitr.resize(maxn);
            for(int i = 1; i < maxn; i *= 2) {
                int ti = i / 2;
                for(int j = 0; j < i; j++) {
                    w[i + j] = std::polar(ftype(1), pi * j / i);
                    if(ti) {
                        bitr[i + j] = 2 * bitr[ti + j % ti] + (j >= ti);
                    }
                }
            }
            initiated = 1;
        }
    }
    
    void fft(auto &a, int n) {
        init();
        if(n == 1) {
            return;
        }
        int hn = n / 2;
        for(int i = 0; i < n; i++) {
            int ti = 2 * bitr[hn + i % hn] + (i > hn);
            if(i < ti) {
                std::swap(a[i], a[ti]);
            }
        }
        for(int i = 1; i < n; i *= 2) {
            for(int j = 0; j < n; j += 2 * i) {
                for(int k = j; k < j + i; k++) {
                    point t = a[k + i] * w[i + k - j];
                    a[k + i] = a[k] - t;
                    a[k] += t;
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
    
        auto operator *= (dft const& B) {
            assert(A.size() == B.A.size());
            size_t n = A.size();
            if(!n) {
                return std::vector<base>();
            }
            for(size_t i = 0; i < n; i++) {
                A[i] *= B[i];
            }
            fft(A, n);
            reverse(begin(A) + 1, end(A));
            std::vector<base> res(n);
            for(size_t i = 0; i < n; i++) {
                res[i] = A[i];
                res[i] /= n;
            }
            return res;
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
    
        auto operator *= (dft const& B) {
            assert(A.size() == B.A.size());
            size_t n = A.size();
            if(!n) {
                return std::vector<base>();
            }
            std::vector<point> C(n);
            for(size_t i = 0; i < n; i++) {
                C[i] = A[i] * (B[i] + conj(B[(n - i) % n]));
                A[i] = A[i] * (B[i] - conj(B[(n - i) % n]));
            }
            fft(C, n);
            fft(A, n);
            reverse(begin(C) + 1, end(C));
            reverse(begin(A) + 1, end(A));
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
