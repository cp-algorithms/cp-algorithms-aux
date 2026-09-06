#ifndef CP_ALGO_MATH_CVECTOR_HPP
#define CP_ALGO_MATH_CVECTOR_HPP
#include "../util/simd.hpp"
#include "../util/complex.hpp"
#include "../util/checkpoint.hpp"
#include "../util/big_alloc.hpp"
#include <ranges>
#include <bit>
CP_ALGO_SIMD_PRAGMA_PUSH
namespace stdx = std::experimental;
namespace cp_algo::math::fft {
    static constexpr size_t flen = 4;
    using ftype = double;
    using vftype = dx4;
    using point = complex<ftype>;
    using vpoint = complex<vftype>;
    static constexpr vftype vz = {};
    vpoint vi(vpoint const& r) {
        return {-imag(r), real(r)};
    }

    struct cvector {
        big_vector<vpoint> r;
        cvector(size_t n) {
            n = std::max(flen, std::bit_ceil(n));
            r.resize(n / flen);
            prepare_roots(n / 16);
            checkpoint("cvector create");
        }

        vpoint& at(size_t k) {return r[k / flen];}
        vpoint at(size_t k) const {return r[k / flen];}
        template<class pt = point>
        inline void set(size_t k, pt const& t) {
            if constexpr(std::is_same_v<pt, point>) {
                real(r[k / flen])[k % flen] = real(t);
                imag(r[k / flen])[k % flen] = imag(t);
            } else {
                at(k) = t;
            }
        }
        template<class pt = point>
        inline pt get(size_t k) const {
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
            } else if(n / 4 - pre_evals < extra.size()) {
                return extra[n / 4 - pre_evals];
            } else {
                return polar<ftype>(1., std::numbers::pi / (ftype)std::bit_floor(n) * (ftype)eval_arg(n));
            }
        }
        static constexpr std::array<point, 32> roots = []() {
            std::array<point, 32> res;
            for(size_t i = 2; i < 32; i++) {
                res[i] = polar<ftype>(1., std::numbers::pi / (1ull << (i - 2)));
            }
            return res;
        }();
        static constexpr point root(size_t n) {
            return roots[std::bit_width(n)];
        }
        template<int step>
        static void exec_on_eval(size_t n, size_t k, auto &&callback) {
            callback(k, root(4 * step * n) * eval_point(step * k));
        }
        template<int step>
        static void exec_on_evals(size_t n, auto &&callback) {
            point factor = root(4 * step * n);
            for(size_t i = 0; i < n; i++) {
                callback(i, factor * eval_point(step * i));
            }
        }

        static void do_dot_iter(point rt, vpoint& Bv, vpoint const& Av, vpoint& res) {
            res += Av * Bv;
            real(Bv) = rotate_right(real(Bv));
            imag(Bv) = rotate_right(imag(Bv));
            auto x = real(Bv)[0], y = imag(Bv)[0];
            real(Bv)[0] = x * real(rt) - y * imag(rt);
            imag(Bv)[0] = x * imag(rt) + y * real(rt);
        }

        void dot(cvector const& t) {
            size_t n = this->size();
            exec_on_evals<1>(n / flen, [&](size_t k, point rt) __attribute__((always_inline)) {
                k *= flen;
                auto [Ax, Ay] = at(k);
                auto Bv = t.at(k);
                vpoint res = vz;
                for (size_t i = 0; i < flen; i++) {
                    vpoint Av = vpoint(vz + Ax[i], vz + Ay[i]);
                    do_dot_iter(rt, Bv, Av, res);
                }
                set(k, res);
            });
            checkpoint("dot");
        }
        // normalize=false leaves the inverse-transform scale for the caller.
        // Only the first count coefficients are promised; the rest are unspecified.
        template<bool partial = true, bool normalize = true>
        void ifft(size_t count = SIZE_MAX) {
            size_t n = size();
            if(count == 0) {return;}
            count = std::min(n, (std::min(count, n) + flen - 1) / flen * flen);
            auto transform = [&]<bool prune>() {
                if constexpr (!partial) {
                    prepare_roots(n / 4);
                    point pi(0, 1);
                    exec_on_evals<4>(n / 4, [&](size_t k, point rt) __attribute__((always_inline)) {
                        k *= 4;
                        point v1 = conj(rt);
                        point v2 = v1 * v1;
                        point v3 = v1 * v2;
                        auto A = get(k);
                        auto B = get(k + 1);
                        auto C = get(k + 2);
                        auto D = get(k + 3);
                        set(k, (A + B) + (C + D));
                        set(k + 2, ((A + B) - (C + D)) * v2);
                        set(k + 1, ((A - B) - pi * (C - D)) * v1);
                        set(k + 3, ((A - B) + pi * (C - D)) * v3);
                    });
                }
                bool parity = std::countr_zero(n) % 2;
                if(parity) {
                    exec_on_evals<2>(n / (2 * flen), [&](size_t k, point rt) __attribute__((always_inline)) {
                        k *= 2 * flen;
                        vpoint cvrt = {vz + real(rt), vz - imag(rt)};
                        auto B = at(k) - at(k + flen);
                        at(k) += at(k + flen);
                        at(k + flen) = B * cvrt;
                    });
                }

                for(size_t leaf = 3 * flen; leaf < n; leaf += 4 * flen) {
                    size_t level = std::countr_one(leaf + 3);
                    for(size_t lvl = 4 + parity; lvl <= level; lvl += 2) {
                        size_t i = (1 << lvl) / 4;
                        exec_on_eval<4>(n >> lvl, leaf >> lvl, [&](size_t k, point rt) __attribute__((always_inline)) {
                            k <<= lvl;
                            vpoint v1 = {vz + real(rt), vz - imag(rt)};
                            vpoint v2 = v1 * v1;
                            vpoint v3 = v1 * v2;
                            for(size_t j = k; j < k + (prune ? std::min(count, i) : i); j += flen) {
                                auto A = at(j);
                                auto B = at(j + i);
                                auto C = at(j + 2 * i);
                                auto D = at(j + 3 * i);
                                at(j) = ((A + B) + (C + D));
                                if(!prune || j - k + 2 * i < count) {at(j + 2 * i) = ((A + B) - (C + D)) * v2;}
                                if(!prune || j - k + i < count) {at(j + i) = ((A - B) - vi(C - D)) * v1;}
                                if(!prune || j - k + 3 * i < count) {at(j + 3 * i) = ((A - B) + vi(C - D)) * v3;}
                            }
                        });
                    }
                }
                checkpoint("ifft");
                if constexpr(normalize) {
                    auto scale = vz + ftype(partial ? flen : 1) / ftype(n);
                    for(size_t k = 0; k < count; k += flen) {
                        set(k, get<vpoint>(k) * scale);
                    }
                }
            };
            // Branching costs more than it saves for nearly full transforms.
            if(count <= n / 4) {transform.template operator()<true>();}
            else {transform.template operator()<false>();}
        }
        // Coefficients at indices >= nonzero must be zero. All evaluations are returned.
        // partial controls the four-lane layout, independently of prefix pruning.
        template<bool partial = true>
        void fft(size_t nonzero = SIZE_MAX) {
            size_t n = size();
            if(nonzero == 0) {return;}
            nonzero = std::min(n, (std::min(nonzero, n) + flen - 1) / flen * flen);
            auto transform = [&]<bool prune>() {
                bool parity = std::countr_zero(n) % 2;
                for(size_t leaf = 0; leaf < n; leaf += 4 * flen) {
                    size_t level = std::countr_zero(n + leaf);
                    level -= level % 2 != parity;
                    for(size_t lvl = level; lvl >= 4; lvl -= 2) {
                        size_t i = (1 << lvl) / 4;
                        exec_on_eval<4>(n >> lvl, leaf >> lvl, [&](size_t k, point rt) __attribute__((always_inline)) {
                            k <<= lvl;
                            vpoint v1 = {vz + real(rt), vz + imag(rt)};
                            vpoint v2 = v1 * v1;
                            vpoint v3 = v1 * v2;
                            auto run = [&]<int terms>(size_t lo, size_t hi) __attribute__((always_inline)) {
                                for(size_t j = k + lo; j < k + hi; j += flen) {
                                    auto A = at(j);
                                    if constexpr(terms == 1) {
                                        at(j + i) = at(j + 2 * i) = at(j + 3 * i) = A;
                                    } else {
                                        auto B = at(j + i) * v1;
                                        auto C = vpoint{}, D = vpoint{};
                                        if constexpr(terms > 2) {C = at(j + 2 * i) * v2;}
                                        if constexpr(terms > 3) {D = at(j + 3 * i) * v3;}
                                        at(j)         = (A + C) + (B + D);
                                        at(j + i)     = (A + C) - (B + D);
                                        at(j + 2 * i) = (A - C) + vi(B - D);
                                        at(j + 3 * i) = (A - C) - vi(B - D);
                                    }
                                }
                            };
                            if constexpr(!prune) {run.template operator()<4>(0, i);}
                            else {
                                size_t q = nonzero >> (lvl - 2), tail = nonzero & (i - 1);
                                if(q >= 4) {run.template operator()<4>(0, i);}
                                else if(q == 3) {run.template operator()<4>(0, tail); run.template operator()<3>(tail, i);}
                                else if(q == 2) {run.template operator()<3>(0, tail); run.template operator()<2>(tail, i);}
                                else if(q == 1) {run.template operator()<2>(0, tail); run.template operator()<1>(tail, i);}
                                else {run.template operator()<1>(0, tail);}
                            }
                        });
                    }
                }
                if(parity) {
                    exec_on_evals<2>(n / (2 * flen), [&](size_t k, point rt) __attribute__((always_inline)) {
                        k *= 2 * flen;
                        vpoint vrt = {vz + real(rt), vz + imag(rt)};
                        auto t = at(k + flen) * vrt;
                        at(k + flen) = at(k) - t;
                        at(k) += t;
                    });
                }
                if constexpr (!partial) {
                    prepare_roots(n / 4);
                    point pi(0, 1);
                    exec_on_evals<4>(n / 4, [&](size_t k, point rt) __attribute__((always_inline)) {
                        k *= 4;
                        point v1 = rt;
                        point v2 = v1 * v1;
                        point v3 = v1 * v2;
                        auto A = get(k);
                        auto B = get(k + 1) * v1;
                        auto C = get(k + 2) * v2;
                        auto D = get(k + 3) * v3;
                        set(k, (A + C) + (B + D));
                        set(k + 1, (A + C) - (B + D));
                        set(k + 2, (A - C) + pi * (B - D));
                        set(k + 3, (A - C) - pi * (B - D));
                    });
                }
                checkpoint("fft");
            };
            if(nonzero <= n / 4) {transform.template operator()<true>();}
            else {transform.template operator()<false>();}
        }
        static constexpr size_t pre_evals = 1 << 16;
        static const std::array<size_t, pre_evals> eval_args;
        static const std::array<point, pre_evals> evalp;
    private:
        static big_vector<point> extra;
        // Keep the usual table small; cache additional roots for large transforms.
        static void prepare_roots(size_t n) {
            if(n <= pre_evals + extra.size()) {return;}
            size_t old = extra.size();
            extra.resize(std::bit_ceil(n) - pre_evals);
            for(size_t i = old; i < extra.size(); i++) {
                size_t j = 4 * (i + pre_evals);
                extra[i] = polar<ftype>(1., std::numbers::pi / (ftype)std::bit_floor(j) * (ftype)eval_arg(j));
            }
        }
    };

    big_vector<point> cvector::extra;

    const std::array<size_t, cvector::pre_evals> cvector::eval_args = []() {
        std::array<size_t, pre_evals> res = {};
        for(size_t i = 1; i < pre_evals; i++) {
            res[i] = res[i >> 1] | (i & 1) << (std::bit_width(i) - 1);
        }
        return res;
    }();
    const std::array<point, cvector::pre_evals> cvector::evalp = []() {
        std::array<point, pre_evals> res = {};
        res[0] = 1;
        for(size_t n = 1; n < pre_evals; n++) {
            res[n] = polar<ftype>(1., std::numbers::pi * ftype(eval_args[n]) / ftype(4 * std::bit_floor(n)));
        }
        return res;
    }();
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_CVECTOR_HPP
