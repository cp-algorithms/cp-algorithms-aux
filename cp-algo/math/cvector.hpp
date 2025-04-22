#ifndef CP_ALGO_MATH_CVECTOR_HPP
#define CP_ALGO_MATH_CVECTOR_HPP
#include "../util/complex.hpp"
#include "../util/checkpoint.hpp"
#include <experimental/simd>
#include <ranges>

namespace stdx = std::experimental;
namespace cp_algo::math::fft {
    using ftype = double;
    static constexpr size_t bytes = 32;
    static constexpr size_t flen = bytes / sizeof(ftype);
    using point = complex<ftype>;
    using vftype [[gnu::vector_size(bytes)]] = ftype;
    using vpoint = complex<vftype>;
    static constexpr vftype fz = {};

    struct cvector {
        std::vector<vftype> x, y;
        cvector(size_t n) {
            n = std::max(flen, std::bit_ceil(n));
            x.resize(n / flen);
            y.resize(n / flen);
            checkpoint("cvector create");
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

        static constexpr size_t pre_roots = 1 << 16;
        static constexpr std::array<point, pre_roots> roots = []() {
            std::array<point, pre_roots> res = {};
            for(size_t n = 1; n < res.size(); n *= 2) {
                for(size_t k = 0; k < n; k++) {
                    res[n + k] = polar(1., std::numbers::pi / ftype(n) * ftype(k));
                }
            }
            return res;
        }();
        static constexpr std::array<size_t, pre_roots> eval_args = []() {
            std::array<size_t, pre_roots> res = {};
            for(size_t i = 1; i < pre_roots; i++) {
                res[i] = res[i >> 1] | (i & 1) << (std::bit_width(i) - 1);
            }
            return res;
        }();
        static constexpr std::array<point, pre_roots> evalp = []() {
            std::array<point, pre_roots> res = {};
            res[0] = 1;
            for(size_t n = 1; n < pre_roots; n++) {
                res[n] = polar(1., std::numbers::pi * ftype(eval_args[n]) / ftype(2 * std::bit_floor(n)));
            }
            return res;
        }();
        static size_t eval_arg(size_t n) {
            if(n < pre_roots) {
                return eval_args[n];
            } else {
                return eval_arg(n / 2) | (n & 1) << (std::bit_width(n) - 1);
            }
        }
        static auto root(size_t n, size_t k) {
            if(n < pre_roots) {
                return roots[n + k];
            } else {
                return polar(1., std::numbers::pi / (ftype)n * (ftype)k);
            }
        }
        static point eval_point(size_t n) {
            if(n < pre_roots) {
                return evalp[n];
            } else {
                return root(2 * std::bit_floor(n), eval_arg(n));
            }
        }
        static void exec_on_roots(size_t n, size_t m, auto &&callback) {
            point cur;
            point arg = root(n, 1);
            for(size_t i = 0; i < m; i++) {
                if(i % 32 == 0 || n < pre_roots) {
                    cur = root(n, i);
                } else {
                    cur *= arg;
                }
                callback(i, cur);
            }
        }
        template<int step = 1>
        static void exec_on_evals(size_t n, auto &&callback) {
            for(size_t i = 0; i < n; i++) {
                callback(i, eval_point(step * i));
            }
        }
        static auto dot_block(size_t k, cvector const& A, cvector const& B) {
            auto rt = eval_point(k / flen / 2);
            if(k / flen % 2) {
                rt = -rt;
            }
            auto [Bvx, Bvy] = B.vget(k);
            auto [Brvx, Brvy] = vpoint(Bvx, Bvy) * vpoint(fz + real(rt), fz + imag(rt));
            auto [Ax, Ay] = A.vget(k);
            vftype Bx[2] = {Brvx, Bvx}, By[2] = {Brvy, Bvy};
            vpoint res = {fz, fz};
            for (size_t i = 0; i < flen; i++) {
                auto Bsx = (vftype*)((ftype*)Bx + flen - i);
                auto Bsy = (vftype*)((ftype*)By + flen - i);
                res += vpoint(fz + Ax[i], fz + Ay[i]) * vpoint{*Bsx, *Bsy};
            }
            return res;
        }

        void dot(cvector const& t) {
            size_t n = this->size();
            for(size_t k = 0; k < n; k += flen) {
                set(k, dot_block(k, *this, t));
            }
            checkpoint("dot");
        }

        void ifft() {
            size_t n = size();
            for(size_t i = flen; i <= n / 2; i *= 2) {
                if (4 * i <= n) { // radix-4
                    exec_on_evals<2>(n / (4 * i), [&](size_t k, point rt) {
                        k *= 4 * i;
                        vpoint v1 = {fz + real(rt), fz - imag(rt)};
                        vpoint v2 = v1 * v1;
                        vpoint v3 = v1 * v2;
                        for(size_t j = k; j < k + i; j += flen) {
                            auto A = get<vpoint>(j);
                            auto B = get<vpoint>(j + i);
                            auto C = get<vpoint>(j + 2 * i);
                            auto D = get<vpoint>(j + 3 * i);
                            set(j        , (A + B + C + D));
                            set(j + 2 * i, (A + B - C - D) * v2);
                            set(j +     i, (A - B - vpoint(fz, fz + 1) * (C - D)) * v1);
                            set(j + 3 * i, (A - B + vpoint(fz, fz + 1) * (C - D)) * v3);
                        }
                    });
                    i *= 2;
                } else { // radix-2 fallback
                    exec_on_evals(n / (2 * i), [&](size_t k, point rt) {
                        k *= 2 * i;
                        vpoint cvrt = {fz + real(rt), fz - imag(rt)};
                        for(size_t j = k; j < k + i; j += flen) {
                            auto A = get<vpoint>(j) + get<vpoint>(j + i);
                            auto B = get<vpoint>(j) - get<vpoint>(j + i);
                            set(j, A);
                            set(j + i, B * cvrt);
                        }
                    });
                }
            }
            checkpoint("ifft");
            for(size_t k = 0; k < n; k += flen) {
                set(k, get<vpoint>(k) /= fz + (ftype)(n / flen));
            }
        }
        void fft() {
            size_t n = size();
            for(size_t i = n / 2; i >= flen; i /= 2) {
                if (i / 2 >= flen) { // radix-4
                    i /= 2;
                    exec_on_evals<2>(n / (4 * i), [&](size_t k, point rt) {
                        k *= 4 * i;
                        vpoint v1 = {fz + real(rt), fz + imag(rt)};
                        vpoint v2 = v1 * v1;
                        vpoint v3 = v1 * v2;
                        for(size_t j = k; j < k + i; j += flen) {
                            auto A = get<vpoint>(j);
                            auto B = get<vpoint>(j + i) * v1;
                            auto C = get<vpoint>(j + 2 * i) * v2;
                            auto D = get<vpoint>(j + 3 * i) * v3;
                            set(j    , (A + C) + (B + D));
                            set(j + i, (A + C) - (B + D));
                            set(j + 2 * i, (A - C) + vpoint(fz, fz + 1) * (B - D));
                            set(j + 3 * i, (A - C) - vpoint(fz, fz + 1) * (B - D));
                        }
                    });
                } else { // radix-2 fallback
                    exec_on_evals(n / (2 * i), [&](size_t k, point rt) {
                        k *= 2 * i;
                        vpoint vrt = {fz + real(rt), fz + imag(rt)};
                        for(size_t j = k; j < k + i; j += flen) {
                            auto t = get<vpoint>(j + i) * vrt;
                            set(j + i, get<vpoint>(j) - t);
                            set(j, get<vpoint>(j) + t);
                        }
                    });
                }
            }
            checkpoint("fft");
        }
    };
}
#endif // CP_ALGO_MATH_CVECTOR_HPP
