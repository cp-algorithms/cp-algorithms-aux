#ifndef CP_ALGO_MATH_CVECTOR_HPP
#define CP_ALGO_MATH_CVECTOR_HPP
#include <algorithm>
#include <complex>
#include <vector>
#include <ranges>
namespace cp_algo::math::fft {
    using ftype = double;
    using point = std::complex<ftype>;

    struct ftvec: std::vector<point> {
        static constexpr size_t pre_roots = 1 << 16;
        static constexpr size_t threshold = 32;
        ftvec(size_t n) {
            this->resize(std::max(threshold, std::bit_ceil(n)));
        }
        static auto dot_block(size_t k, ftvec const& A, ftvec const& B) {
            static std::array<point, 2 * threshold> r;
            std::ranges::fill(r, point(0));
            for(size_t i = 0; i < threshold; i++) {
                for(size_t j = 0; j < threshold; j++) {
                    r[i + j] += A[k + i] * B[k + j];
                }
            }
            auto rt = ftype(k / threshold % 2 ? -1 : 1) * eval_point(k / threshold / 2);
            static std::array<point, threshold> res;
            for(size_t i = 0; i < threshold; i++) {
                res[i] = r[i] + r[i + threshold] * rt;
            }
            return res;
        }

        void dot(ftvec const& t) {
            size_t n = this->size();
            for(size_t k = 0; k < n; k += threshold) {
                std::ranges::copy(dot_block(k, *this, t), this->begin() + k);
            }
        }
        static std::array<point, pre_roots> roots, evalp;
        static std::array<size_t, pre_roots> eval_args;
        static point root(size_t n, size_t k) {
            if(n + k < pre_roots && roots[n + k] != point{}) {
                return roots[n + k];
            }
            auto res = std::polar(1., std::numbers::pi * ftype(k) / ftype(n));
            if(n + k < pre_roots) {
                roots[n + k] = res;
            }
            return res;
        }
        static size_t eval_arg(size_t n) {
            if(n < pre_roots && eval_args[n]) {
                return eval_args[n];
            } else if(n == 0) {
                return 0;
            }
            auto res = eval_arg(n / 2) | (n & 1) << (std::bit_width(n) - 1);
            if(n < pre_roots) {
                eval_args[n] = res;
            }
            return res;
        }
        static point eval_point(size_t n) {
            if(n < pre_roots && evalp[n] != point{}) {
                return evalp[n];
            } else if(n == 0) {
                return point(1);
            }
            auto res = root(2 * std::bit_floor(n), eval_arg(n));
            if(n < pre_roots) {
                evalp[n] = res;
            }
            return res;
        }
        static void exec_on_roots(size_t n, size_t m, auto &&callback) {
            auto step = root(n, 1);
            auto rt = point(1);
            for(size_t i = 0; i < m; i++) {
                if(i % threshold == 0) {
                    rt = root(n / threshold, i / threshold);
                }
                callback(i, rt);
                rt *= step;
            }
        }
        static void exec_on_evals(size_t n, auto &&callback) {
            for(size_t i = 0; i < n; i++) {
                callback(i, eval_point(i));
            }
        }

        void ifft() {
            size_t n = this->size();
            for(size_t half = threshold; half <= n / 2; half *= 2) {
                exec_on_evals(n / (2 * half), [&](size_t k, point rt) {
                    k *= 2 * half;
                    for(size_t j = k; j < k + half; j++) {
                        auto A = this->at(j) + this->at(j + half);
                        auto B = this->at(j) - this->at(j + half);
                        this->at(j) = A;
                        this->at(j + half) = B * conj(rt);
                    }
                });
            }
            point ni = point(int(threshold)) / point(int(n));
            for(auto &it: *this) {
                it *= ni;
            }
        }
        void fft() {
            size_t n = this->size();
            for(size_t half = n / 2; half >= threshold; half /= 2) {
                exec_on_evals(n / (2 * half), [&](size_t k, point rt) {
                    k *= 2 * half;
                    for(size_t j = k; j < k + half; j++) {
                        auto t = this->at(j + half) * rt;
                        this->at(j + half) = this->at(j) - t;
                        this->at(j) += t;
                    }
                });
            }
        }
    };
    std::array<point, ftvec::pre_roots> ftvec::roots = {};
    std::array<point, ftvec::pre_roots> ftvec::evalp = {};
    std::array<size_t, ftvec::pre_roots> ftvec::eval_args = {};
}
#endif // CP_ALGO_MATH_CVECTOR_HPP
