#ifndef CP_ALGO_MATH_CVECTOR_HPP
#define CP_ALGO_MATH_CVECTOR_HPP
#include "../util/simd.hpp"
#include "../util/complex.hpp"
#include "../util/checkpoint.hpp"
#include "../util/big_alloc.hpp"
#include <ranges>

namespace stdx = std::experimental;
namespace cp_algo::math::fft {
    static constexpr size_t flen = 4;
    using ftype = double;
    using vftype = simd<ftype, flen>;
    using point = complex<ftype>;
    using vpoint = complex<vftype>;
    static constexpr vftype vz = {};
    vpoint vi(vpoint const& r) {
        return {-imag(r), real(r)};
    }

    struct cvector {
        std::vector<vpoint, big_alloc<vpoint>> r;
        cvector(size_t n) {
            n = std::max(flen, std::bit_ceil(n));
            r.resize(n / flen);
            checkpoint("cvector create");
        }

        vpoint& at(size_t k) {return r[k / flen];}
        vpoint at(size_t k) const {return r[k / flen];}
        template<class pt = point>
        void set(size_t k, pt t) {
            if constexpr(std::is_same_v<pt, point>) {
                real(r[k / flen])[k % flen] = real(t);
                imag(r[k / flen])[k % flen] = imag(t);
            } else {
                at(k) = t;
            }
        }
        template<class pt = point>
        pt get(size_t k) const {
            if constexpr(std::is_same_v<pt, point>) {
                return {real(r[k / flen])[k % flen], imag(r[k / flen])[k % flen]};
            } else {
                return at(k);
            }
        }

        size_t size() const {
            return flen * r.size();
        }
        static constexpr size_t eval_arg(size_t n) {
            if(n < pre_evals) {
                return eval_args[n];
            } else {
                return eval_arg(n / 2) | (n & 1) << (std::bit_width(n) - 1);
            }
        }
        static constexpr point eval_point(size_t n) {
            if(n % 2) {
                return -eval_point(n - 1);
            } else if(n % 4) {
                return eval_point(n - 2) * point(0, 1);
            } else if(n / 4 < pre_evals) {
                return evalp[n / 4];
            } else {
                return polar(1., std::numbers::pi / (ftype)std::bit_floor(n) * (ftype)eval_arg(n));
            }
        }
        static constexpr point root(size_t n) {
            return eval_point(n / 2);
        }
        template<int step>
        static void exec_on_evals(size_t n, auto &&callback) {
            point factor = root(4 * step * n);
            for(size_t i = 0; i < n; i++) {
                callback(i, factor * eval_point(step * i));
            }
        }
        template<int step>
        static void exec_on_eval(size_t n, size_t k, auto &&callback) {
            callback(root(4 * step * n) * eval_point(step * k));
        }

        void dot(cvector const& t) {
            size_t n = this->size();
            exec_on_evals<1>(n / flen, [&](size_t k, point rt) {
                k *= flen;
                auto [Ax, Ay] = at(k);
                auto Bv = t.at(k);
                vpoint res = vz;
                for (size_t i = 0; i < flen; i++) {
                    res += vpoint(vz + Ax[i], vz + Ay[i]) * Bv;
                    real(Bv) = __builtin_shufflevector(real(Bv), real(Bv), 3, 0, 1, 2);
                    imag(Bv) = __builtin_shufflevector(imag(Bv), imag(Bv), 3, 0, 1, 2);
                    auto x = real(Bv)[0], y = imag(Bv)[0];
                    real(Bv)[0] = x * real(rt) - y * imag(rt);
                    imag(Bv)[0] = x * imag(rt) + y * real(rt);
                }
                set(k, res);
            });
            checkpoint("dot");
        }

        void ifft() {
            size_t n = size();
            for(size_t i = flen; i <= n / 2; i *= 2) {
                if (4 * i <= n) { // radix-4
                    exec_on_evals<4>(n / (4 * i), [&](size_t k, point rt) {
                        k *= 4 * i;
                        vpoint v1 = {vz + real(rt), vz - imag(rt)};
                        vpoint v2 = v1 * v1;
                        vpoint v3 = v1 * v2;
                        for(size_t j = k; j < k + i; j += flen) {
                            auto A = at(j);
                            auto B = at(j + i);
                            auto C = at(j + 2 * i);
                            auto D = at(j + 3 * i);
                            at(j) = (A + B + C + D);
                            at(j + 2 * i) = (A + B - C - D) * v2;
                            at(j +     i) = (A - B - vi(C - D)) * v1;
                            at(j + 3 * i) = (A - B + vi(C - D)) * v3;
                        }
                    });
                    i *= 2;
                } else { // radix-2 fallback
                    exec_on_evals<2>(n / (2 * i), [&](size_t k, point rt) {
                        k *= 2 * i;
                        vpoint cvrt = {vz + real(rt), vz - imag(rt)};
                        for(size_t j = k; j < k + i; j += flen) {
                            auto B = at(j) - at(j + i);
                            at(j) += at(j + i);
                            at(j + i) = B * cvrt;
                        }
                    });
                }
            }
            checkpoint("ifft");
            for(size_t k = 0; k < n; k += flen) {
                set(k, get<vpoint>(k) /= vz + (ftype)(n / flen));
            }
        }
        void fft() {
            size_t n = size();
            for(size_t i = n / 2; i >= flen; i /= 2) {
                if (i / 2 >= flen) { // radix-4
                    i /= 2;
                    exec_on_evals<4>(n / (4 * i), [&](size_t k, point rt) {
                        k *= 4 * i;
                        vpoint v1 = {vz + real(rt), vz + imag(rt)};
                        vpoint v2 = v1 * v1;
                        vpoint v3 = v1 * v2;
                        for(size_t j = k; j < k + i; j += flen) {
                            auto A = at(j);
                            auto B = at(j + i) * v1;
                            auto C = at(j + 2 * i) * v2;
                            auto D = at(j + 3 * i) * v3;
                            at(j)         = (A + C) + (B + D);
                            at(j + i)     = (A + C) - (B + D);
                            at(j + 2 * i) = (A - C) + vi(B - D);
                            at(j + 3 * i) = (A - C) - vi(B - D);
                        }
                    });
                } else { // radix-2 fallback
                    exec_on_evals<2>(n / (2 * i), [&](size_t k, point rt) {
                        k *= 2 * i;
                        vpoint vrt = {vz + real(rt), vz + imag(rt)};
                        for(size_t j = k; j < k + i; j += flen) {
                            auto t = at(j + i) * vrt;
                            at(j + i) = at(j) - t;
                            at(j) += t;
                        }
                    });
                }
            }
            checkpoint("fft");
        }
        static constexpr size_t pre_evals = 1 << 16;
        static constexpr std::array<size_t, pre_evals> eval_args = []() {
            std::array<size_t, pre_evals> res = {};
            for(size_t i = 1; i < pre_evals; i++) {
                res[i] = res[i >> 1] | (i & 1) << (std::bit_width(i) - 1);
            }
            return res;
        }();
        static constexpr std::array<point, pre_evals> evalp = []() {
            std::array<point, pre_evals> res = {};
            res[0] = 1;
            for(size_t n = 1; n < pre_evals; n++) {
                res[n] = polar(1., std::numbers::pi * ftype(eval_args[n]) / ftype(4 * std::bit_floor(n)));
            }
            return res;
        }();
    };
}
#endif // CP_ALGO_MATH_CVECTOR_HPP
