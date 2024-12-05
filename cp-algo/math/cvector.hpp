#ifndef CP_ALGO_MATH_CVECTOR_HPP
#define CP_ALGO_MATH_CVECTOR_HPP
#include <algorithm>
#include <cassert>
#include <complex>
#include <vector>
#include <ranges>
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
                auto arg = std::numbers::pi / ftype(n);
                if constexpr(std::is_same_v<pt, point>) {
                    return {cos(ftype(k) * arg), sin(ftype(k) * arg)};
                } else {
                    return WITH_IV(pt{vftype{cos(ftype(k + i) * arg)...},
                                      vftype{sin(ftype(k + i) * arg)...}});
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
                set(k, get<vpoint>(k) /= to_pt<vpoint>(ftype(n)));
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
            auto base = std::polar(1., std::numbers::pi / ftype(n));
            point cur = 1;
            for(size_t k = 0; k < n; k++) {
                if((k & 15) == 0) {
                    cur = std::polar(1., std::numbers::pi * ftype(k) / ftype(n));
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
}
#endif // CP_ALGO_MATH_CVECTOR_HPP
