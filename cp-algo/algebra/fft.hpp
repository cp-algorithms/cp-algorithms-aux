#ifndef ALGEBRA_FFT_HPP
#define ALGEBRA_FFT_HPP
#include "common.hpp"
#include "modular.hpp"
#include <algorithm>
#include "cassert"
#include <vector>
namespace algebra {
    namespace fft {
        using ftype = double;
        struct point {
            ftype x, y;
            
            ftype real() {return x;}
            ftype imag() {return y;}
            
            point(): x(0), y(0){}
            point(ftype x, ftype y = 0): x(x), y(y){}
            
            static point polar(ftype rho, ftype ang) {
                return point{rho * cos(ang), rho * sin(ang)};
            }
            
            point conj() const {
                return {x, -y};
            }
            
            point operator +=(const point &t) {x += t.x, y += t.y; return *this;}
            point operator +(const point &t) const {return point(*this) += t;}
            point operator -(const point &t) const {return {x - t.x, y - t.y};}
            point operator *(const point &t) const {return {x * t.x - y * t.y, x * t.y + y * t.x};}
        };

        point w[maxn]; // w[2^n + k] = exp(pi * k / (2^n))
        int bitr[maxn];// b[2^n + k] = bitreverse(k)
        const ftype pi = acos(-1);
        bool initiated = 0;
        void init() {
            if(!initiated) {
                for(int i = 1; i < maxn; i *= 2) {
                    int ti = i / 2;
                    for(int j = 0; j < i; j++) {
                        w[i + j] = point::polar(ftype(1), pi * j / i);
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
        
        void mul_slow(std::vector<auto> &a, const std::vector<auto> &b) {
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
        
        template<int m>
        struct dft {
            static constexpr int split = 1 << 15;
            std::vector<point> A;
            
            dft(std::vector<modular<m>> const& a, size_t n): A(n) {
                for(size_t i = 0; i < std::min(n, a.size()); i++) {
                    A[i] = point(
                        a[i].rem() % split,
                        a[i].rem() / split
                    );
                }
                if(n) {
                    fft(A, n);
                }
            }
        
            auto operator * (dft const& B) {
                assert(A.size() == B.A.size());
                size_t n = A.size();
                if(!n) {
                    return std::vector<modular<m>>();
                }
                std::vector<point> C(n), D(n);
                for(size_t i = 0; i < n; i++) {
                    C[i] = A[i] * (B[i] + B[(n - i) % n].conj());
                    D[i] = A[i] * (B[i] - B[(n - i) % n].conj());
                }
                fft(C, n);
                fft(D, n);
                reverse(begin(C) + 1, end(C));
                reverse(begin(D) + 1, end(D));
                int t = 2 * n;
                std::vector<modular<m>> res(n);
                for(size_t i = 0; i < n; i++) {
                    modular<m> A0 = llround(C[i].real() / t);
                    modular<m> A1 = llround(C[i].imag() / t + D[i].imag() / t);
                    modular<m> A2 = llround(D[i].real() / t);
                    res[i] = A0 + A1 * split - A2 * split * split;
                }
                return res;
            }
            
            point& operator [](int i) {return A[i];}
            point operator [](int i) const {return A[i];}
        };
        
        size_t com_size(size_t as, size_t bs) {
            if(!as || !bs) {
                return 0;
            }
            size_t n = as + bs - 1;
            while(__builtin_popcount(n) != 1) {
                n++;
            }
            return n;
        }
        
        template<int m>
        void mul(std::vector<modular<m>> &a, std::vector<modular<m>> b) {
            if(std::min(a.size(), b.size()) < magic) {
                mul_slow(a, b);
                return;
            }
            auto n = com_size(a.size(), b.size());
            auto A = dft<m>(a, n);
            if(a == b) {
                a = A * A;
            } else {
                a = A * dft<m>(b, n);
            }
        }
    }
}
#endif // ALGEBRA_FFT_HPP
