#ifndef CP_ALGO_MATH_FFT_HPP
#define CP_ALGO_MATH_FFT_HPP
#include "common.hpp"
#include "cp-algo/number_theory/modint.hpp"
#include <algorithm>
#include <complex>
#include <cassert>
#include <ranges>
#include <vector>
#include <bit>

namespace cp_algo::math::fft {
    using ftype = double;
    static constexpr size_t bytes = 32;
    static constexpr size_t flen = bytes / sizeof(ftype);
    using point = std::complex<ftype>;
    using vftype [[gnu::vector_size(bytes)]] = ftype;
    using vpoint = std::complex<vftype>;

#define WITH_IV(...)                             \
  [&]<size_t ... i>(std::index_sequence<i...>) { \
      return __VA_ARGS__;                        \
  }(std::make_index_sequence<flen>());

    template<typename ft>
    constexpr ft to_ft(auto x) {
        return ft{} + x;
    }
    template<typename pt>
    constexpr pt to_pt(point r) {
        using ft = std::conditional_t<std::is_same_v<point, pt>, ftype, vftype>;
        return {to_ft<ft>(r.real()), to_ft<ft>(r.imag())};
    }
    struct cvector {
        static constexpr size_t pre_roots = 1 << 17;
        std::vector<vftype> x, y;
        cvector(size_t n) {
            n = std::max(flen, std::bit_ceil(n));
            x.resize(n / flen);
            y.resize(n / flen);
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
        vpoint vget(size_t k) const {
            return get<vpoint>(k);
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
        template<class pt = point>
        static pt root(size_t n, size_t k) {
            if(n < pre_roots) {
                return roots.get<pt>(n + k);
            } else {
                auto arg = std::numbers::pi / n;
                if constexpr(std::is_same_v<pt, point>) {
                    return {cos(k * arg), sin(k * arg)};
                } else {
                    return WITH_IV(pt{vftype{cos((k + i) * arg)...},
                                      vftype{sin((k + i) * arg)...}});
                }
            }
        }
        template<class pt = point>
        static void exec_on_roots(size_t n, size_t m, auto &&callback) {
            size_t step = sizeof(pt) / sizeof(point);
            pt cur;
            pt arg = to_pt<pt>(root<point>(n, step));
            for(size_t i = 0; i < m; i += step) {
                if(i % 64 == 0 || n < pre_roots) {
                    cur = root<pt>(n, i);
                } else {
                    cur *= arg;
                }
                callback(i, cur);
            }
        }

        void ifft() {
            size_t n = size();
            for(size_t i = 1; i < n; i *= 2) {
                for(size_t j = 0; j < n; j += 2 * i) {
                    auto butterfly = [&]<class pt>(size_t k, pt rt) {
                        k += j;
                        auto t = get<pt>(k + i) * conj(rt);
                        set(k + i, get<pt>(k) - t);
                        set(k, get<pt>(k) + t);
                    };
                    if(2 * i <= flen) {
                        exec_on_roots(i, i, butterfly);
                    } else {
                        exec_on_roots<vpoint>(i, i, butterfly);
                    }
                }
            }
            for(size_t k = 0; k < n; k += flen) {
                set(k, get<vpoint>(k) /= to_pt<vpoint>(n));
            }
        }
        void fft() {
            size_t n = size();
            for(size_t i = n / 2; i >= 1; i /= 2) {
                for(size_t j = 0; j < n; j += 2 * i) {
                    auto butterfly = [&]<class pt>(size_t k, pt rt) {
                        k += j;
                        auto A = get<pt>(k) + get<pt>(k + i);
                        auto B = get<pt>(k) - get<pt>(k + i);
                        set(k, A);
                        set(k + i, B * rt);
                    };
                    if(2 * i <= flen) {
                        exec_on_roots(i, i, butterfly);
                    } else {
                        exec_on_roots<vpoint>(i, i, butterfly);
                    }
                }
            }
        }
    };
    const cvector cvector::roots = []() {
        cvector res(pre_roots);
        for(size_t n = 1; n < res.size(); n *= 2) {
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
        int split;
        cvector A, B;
        
        dft(auto const& a, size_t n): A(n), B(n) {
            split = std::sqrt(base::mod());
            cvector::exec_on_roots(2 * n, size(a), [&](size_t i, point rt) {
                size_t ti = std::min(i, i - n);
                A.set(ti, A.get(ti) + ftype(a[i].rem() % split) * rt);
                B.set(ti, B.get(ti) + ftype(a[i].rem() / split) * rt);
    
            });
            if(n) {
                A.fft();
                B.fft();
            }
        }

        void mul(auto &&C, auto const& D, auto &res, size_t k) {
            assert(A.size() == C.size());
            size_t n = A.size();
            if(!n) {
                res = {};
                return;
            }
            for(size_t i = 0; i < n; i += flen) {
                auto tmp = A.vget(i) * D.vget(i) + B.vget(i) * C.vget(i);
                A.set(i, A.vget(i) * C.vget(i));
                B.set(i, B.vget(i) * D.vget(i));
                C.set(i, tmp);
            }
            A.ifft();
            B.ifft();
            C.ifft();
            auto splitsplit = (base(split) * split).rem();
            cvector::exec_on_roots(2 * n, std::min(n, k), [&](size_t i, point rt) {
                rt = conj(rt);
                auto Ai = A.get(i) * rt;
                auto Bi = B.get(i) * rt;
                auto Ci = C.get(i) * rt;
                int64_t A0 = llround(real(Ai));
                int64_t A1 = llround(real(Ci));
                int64_t A2 = llround(real(Bi));
                res[i] = A0 + A1 * split + A2 * splitsplit;
                if(n + i >= k) {
                    return;
                }
                int64_t B0 = llround(imag(Ai));
                int64_t B1 = llround(imag(Ci));
                int64_t B2 = llround(imag(Bi));
                res[n + i] = B0 + B1 * split + B2 * splitsplit;
            });
        }
        void mul_inplace(auto &&B, auto& res, size_t k) {
            mul(B.A, B.B, res, k);
        }
        void mul(auto const& B, auto& res, size_t k) {
            mul(cvector(B.A), B.B, res, k);
        }
        std::vector<base> operator *= (dft &B) {
            std::vector<base> res(2 * A.size());
            mul_inplace(B, res, size(res));
            return res;
        }
        std::vector<base> operator *= (dft const& B) {
            std::vector<base> res(2 * A.size());
            mul(B, res, size(res));
            return res;
        }
        auto operator * (dft const& B) const {
            return dft(*this) *= B;
        }
        
        point operator [](int i) const {return A.get(i);}
    };
    
    void mul_slow(auto &a, auto const& b, size_t k) {
        if(empty(a) || empty(b)) {
            a.clear();
        } else {
            int n = std::min(k, size(a));
            int m = std::min(k, size(b));
            a.resize(k);
            for(int j = k - 1; j >= 0; j--) {
                a[j] *= b[0];
                for(int i = std::max(j - n, 0) + 1; i < std::min(j + 1, m); i++) {
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
        if(std::min({k, size(a), size(b)}) < 64) {
            mul_slow(a, b, k);
            return;
        }
        auto n = std::max(flen, std::bit_ceil(
            std::min(k, size(a)) + std::min(k, size(b)) - 1
        ) / 2);
        a.resize(k);
        auto A = dft<base>(a, n);
        if(&a == &b) {
            A.mul(A, a, k);
        } else {
            A.mul_inplace(dft<base>(std::views::take(b, k), n), a, k);
        }
    }
    void mul(auto &a, auto const& b) {
        if(size(a)) {
            mul_truncate(a, b, size(a) + size(b) - 1);
        }
    }
}
#endif // CP_ALGO_MATH_FFT_HPP
