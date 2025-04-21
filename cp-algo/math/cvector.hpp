#ifndef CP_ALGO_MATH_CVECTOR_HPP
#define CP_ALGO_MATH_CVECTOR_HPP
#include "../util/complex.hpp"
#include "../util/checkpoint.hpp"
#include <experimental/simd>
#include <ranges>

namespace stdx = std::experimental;
namespace cp_algo::math::fft {
    using ftype = double;
    using point = complex<ftype>;
    using vftype = stdx::native_simd<ftype>;
    using vpoint = complex<vftype>;
    static constexpr size_t flen = vftype::size();

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
        static void exec_on_evals(size_t n, auto &&callback) {
            for(size_t i = 0; i < n; i++) {
                callback(i, eval_point(i));
            }
        }
        static auto dot_block(size_t k, cvector const& A, cvector const& B) {
            auto rt = eval_point(k / flen / 2);
            if(k / flen % 2) {
                rt = -rt;
            }
            auto [Bvx, Bvy] = B.vget(k);
            auto [Brvx, Brvy] = vpoint(Bvx, Bvy) * vpoint(real(rt), imag(rt));
            auto [Ax, Ay] = A.vget(k);
            alignas(32) ftype Bx[2 * flen];
            alignas(32) ftype By[2 * flen];
            Bvx.copy_to(Bx + flen, stdx::vector_aligned);
            Bvy.copy_to(By + flen, stdx::vector_aligned);
            Brvx.copy_to(Bx, stdx::vector_aligned);
            Brvy.copy_to(By, stdx::vector_aligned);
            vpoint res = {0, 0};
            for(size_t i = 0; i < flen; i++) {
                vftype Bsx, Bsy;
                Bsx.copy_from(Bx + flen - i, stdx::element_aligned);
                Bsy.copy_from(By + flen - i, stdx::element_aligned);
                res += vpoint(Ax[i], Ay[i]) * vpoint(Bsx, Bsy);
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
                exec_on_evals(n / (2 * i), [&](size_t k, point rt) {
                    k *= 2 * i;
                    vpoint vrt = {real(rt), imag(rt)};
                    for(size_t j = k; j < k + i; j += flen) {
                        auto A = get<vpoint>(j) + get<vpoint>(j + i);
                        auto B = get<vpoint>(j) - get<vpoint>(j + i);
                        set(j, A);
                        set(j + i, B * conj(vrt));
                    }
                });
            }
            checkpoint("ifft");
            for(size_t k = 0; k < n; k += flen) {
                set(k, get<vpoint>(k) /= (ftype)(n / flen));
            }
        }
        void fft() {
            size_t n = size();
            for(size_t i = n / 2; i >= flen; i /= 2) {
                exec_on_evals(n / (2 * i), [&](size_t k, point rt) {
                    k *= 2 * i;
                    vpoint vrt = {real(rt), imag(rt)};
                    for(size_t j = k; j < k + i; j += flen) {
                        auto t = get<vpoint>(j + i) * vrt;
                        set(j + i, get<vpoint>(j) - t);
                        set(j, get<vpoint>(j) + t);
                    }
                });
            }
            checkpoint("fft");
        }
    };
}
#endif // CP_ALGO_MATH_CVECTOR_HPP
