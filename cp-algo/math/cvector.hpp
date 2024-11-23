#ifndef CP_ALGO_MATH_CVECTOR_HPP
#define CP_ALGO_MATH_CVECTOR_HPP
#include "../util/complex.hpp"
#include <experimental/simd>
namespace cp_algo::math::fft {
    using ftype = double;
    using point = complex<ftype>;
    using vftype = std::experimental::native_simd<ftype>;
    using vpoint = complex<vftype>;
    static constexpr size_t flen = vftype::size();

    struct cvector {
        static constexpr size_t pre_roots = 1 << 18;
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
        template< bool precalc = false, class ft = point>
        static auto root(size_t n, size_t k, ft &&arg) {
            if(n < pre_roots && !precalc) {
                return roots.get<complex<ft>>(n + k);
            } else {
                return complex<ft>::polar(1., arg);
            }
        }
        template<class pt = point, bool precalc = false>
        static void exec_on_roots(size_t n, size_t m, auto &&callback) {
            ftype arg = std::numbers::pi / (ftype)n;
            size_t step = sizeof(pt) / sizeof(point);
            using ft = pt::value_type;
            auto k = [&]() {
                if constexpr(std::is_same_v<pt, point>) {
                    return ft{};
                } else {
                    return ft{[](auto i) {return i;}};
                }
            }();
            for(size_t i = 0; i < m; i += step, k += (ftype)step) {
                callback(i, root<precalc>(n, i, arg * k));
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
                    if(i < flen) {
                        exec_on_roots<point>(i, i, butterfly);
                    } else {
                        exec_on_roots<vpoint>(i, i, butterfly);
                    }
                }
            }
            for(size_t k = 0; k < n; k += flen) {
                set(k, get<vpoint>(k) /= (ftype)n);
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
                    if(i < flen) {
                        exec_on_roots<point>(i, i, butterfly);
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
            auto propagate = [&](size_t k, auto rt) {
                res.set(n + k, rt);
            };
            if(n < flen) {
                res.exec_on_roots<point, true>(n, n, propagate);
            } else {
                res.exec_on_roots<vpoint, true>(n, n, propagate);
            }
        }
        return res;
    }();
}
#endif // CP_ALGO_MATH_CVECTOR_HPP
