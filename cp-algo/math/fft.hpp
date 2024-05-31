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
    const auto bitr = [](){
        std::vector<size_t> bitr(maxn);
        for(size_t n = 2; n < maxn; n *= 2) {
            for(size_t k = 0; k < n; k++) {
                bitr[n + k] = bitr[n + k / 2] / 2 + (k & 1) * (n / 2);
            }
        }
        return bitr;
    }();
    size_t bitreverse(size_t n, size_t k) {
        size_t hn = n / 2;
        if(k >= hn && n > 1) {
            return 2 * bitr[k] + 1;
        } else {
            return 2 * bitr[hn + k];
        }
    }

    using ftype = double;
    static constexpr size_t bytes = 32;
    static constexpr size_t flen = bytes / sizeof(ftype);
    using point = std::complex<ftype>;
    using vftype [[gnu::vector_size(bytes)]] = ftype;
    using vpoint = std::complex<vftype>;

    constexpr vftype to_vec(ftype x) {
        return vftype{} + x;
    }
    constexpr vpoint to_vec(point r) {
        return {to_vec(r.real()), to_vec(r.imag())};
    }
    struct cvector {
        std::vector<vftype> x, y;
        cvector() {}
        cvector(size_t n) {
            resize(n);
        }
        void resize(size_t n) {
            n = std::bit_ceil(std::max<size_t>(n, 4));
            if(size() != n) {
                x.resize(n / flen);
                y.resize(n / flen);
            }
        }
        template<class pt = point>
        void set(size_t k, pt t) {
            if constexpr(std::is_same_v<pt, point>) {
                x[k / flen][k % flen] = real(t);
                y[k / flen][k % flen] = imag(t);
            } else {
                x[k / flen] = real(t);
                y[k / flen] = imag(t);
            }
        }
        template<class pt = point>
        pt get(size_t k) const {
            if constexpr(std::is_same_v<pt, point>) {
                return {x[k / flen][k % flen], y[k / flen][k % flen]};
            } else {
                return {x[k / flen], y[k / flen]};
            }
        }
        size_t size() const {
            return flen * std::size(x);
        }
        void dot(cvector const& t) {
            size_t n = size();
            for(size_t k = 0; k < n; k += flen) {
                set(k, get<vpoint>(k) * t.get<vpoint>(k));
            }
        }
        static const cvector roots;

        void ifft() {
            size_t n = size();
            for(size_t i = 1; i < n; i *= 2) {
                for(size_t j = 0; j < n; j += 2 * i) {
                    auto butterfly = [&]<class pt>(pt) {
                        size_t step = sizeof(pt) / sizeof(point);
                        for(size_t k = j; k < j + i; k += step) {
                            auto T = get<pt>(k + i) * conj(roots.get<pt>(i + k - j));
                            set(k + i, get<pt>(k) - T);
                            set(k, get<pt>(k) + T);
                        }
                    };
                    if(2 * i <= flen) {
                        butterfly(point{});
                    } else {
                        butterfly(vpoint{});
                    }
                }
            }
            for(size_t k = 0; k < n; k += flen) {
                set(k, get<vpoint>(k) /= to_vec(n));
            }
        }
        void fft() {
            size_t n = size();
            for(size_t i = n / 2; i >= 1; i /= 2) {
                for(size_t j = 0; j < n; j += 2 * i) {
                    auto butterfly = [&]<class pt>(pt) {
                        size_t step = sizeof(pt) / sizeof(point);
                        for(size_t k = j; k < j + i; k += step) {
                            auto A = get<pt>(k) + get<pt>(k + i);
                            auto B = get<pt>(k) - get<pt>(k + i);
                            set(k, A);
                            set(k + i, B * roots.get<pt>(i + k - j));
                        }
                    };
                    if(2 * i <= flen) {
                        butterfly(point{});
                    } else {
                        butterfly(vpoint{});
                    }
                }
            }
        }
    };
    const cvector cvector::roots = []() {
        cvector res(maxn);
        for(size_t n = 1; n < maxn; n *= 2) {
            auto base = std::polar(1., std::numbers::pi / n);
            point cur = 1;
            for(size_t k = 0; k < n; k++) {
                if((k & 15) == 0) {
                    cur = std::polar(1., std::numbers::pi * k / n);
                }
                res.set(n + k, cur);
                cur *= base;
            }
        }
        return res;
    }();

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
        cvector A;
        
        dft(std::vector<base> const& a, size_t n): A(n) {
            for(size_t i = 0; i < std::min(n, a.size()); i++) {
                A.set(i, a[i]);
            }
            if(n) {
                A.fft();
            }
        }
    
        std::vector<base> operator *= (dft const& B) {
            assert(A.size() == B.A.size());
            size_t n = A.size();
            if(!n) {
                return std::vector<base>();
            }
            A.dot(B.A);
            A.ifft();
            std::vector<base> res(n);
            for(size_t k = 0; k < n; k++) {
                res[k] = A.get(k);
            }
            return res;
        }

        auto operator * (dft const& B) const {
            return dft(*this) *= B;
        }
        
        point operator [](int i) const {return A.get(i);}
    };

    template<modint_type base>
    struct dft<base> {
        static constexpr int split = 1 << 15;
        cvector A;
        
        dft(std::vector<base> const& a, size_t n): A(n) {
            for(size_t i = 0; i < std::min(n, a.size()); i++) {
                A.set(i, point{
                    ftype(a[i].rem() % split),
                    ftype(a[i].rem() / split)
                });
            }
            if(n) {
                A.fft();
            }
        }

        std::vector<base> mul(auto &&B) {
            assert(A.size() == B.size());
            size_t n = A.size();
            if(!n) {
                return std::vector<base>();
            }
            for(size_t i = 0; 2 * i <= n; i++) {
                size_t j = (n - i) % n;
                size_t x = bitreverse(n, i);
                size_t y = bitreverse(n, j);
                auto Ax = A.get(x), Bx = B.get(x);
                auto Ay = A.get(y), By = B.get(y);
                B.set(x, Ax * (Bx + conj(By)));
                A.set(x, Ax * (Bx - conj(By)));
                B.set(y, Ay * (By + conj(Bx)));
                A.set(y, Ay * (By - conj(Bx)));
            }
            A.ifft();
            B.ifft();
            std::vector<base> res(n);
            for(size_t i = 0; i < n; i++) {
                base A0 = llround(B.get(i).real()) / 2;
                base A1 = llround(B.get(i).imag() + A.get(i).imag()) / 2;
                base A2 = llround(A.get(i).real()) / 2;
                res[i] = A0 + A1 * split - A2 * split * split;
            }
            return res;
        }
        std::vector<base> operator *= (auto &&B) {
            return mul(B.A);
        }

        auto operator * (dft const& B) const {
            return dft(*this) *= dft(B);
        }
        
        point operator [](int i) const {return A.get(i);}
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
            a = A *= dft<base>(A);
        } else {
            a = A *= dft<base>(b, n);
        }
    }
    template<typename base>
    void circular_mul(std::vector<base> &a, std::vector<base> const& b) {
        auto n = std::bit_ceil(a.size());
        auto A = dft<base>(a, n);
        if(a == b) {
            a = A *= dft<base>(A);
        } else {
            a = A *= dft<base>(b, n);
        }
    }
}
#endif // CP_ALGO_MATH_FFT_HPP
