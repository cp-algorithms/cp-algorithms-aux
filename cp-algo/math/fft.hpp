#ifndef CP_ALGO_MATH_FFT_HPP
#define CP_ALGO_MATH_FFT_HPP
#include "../number_theory/modint.hpp"
#include "../util/checkpoint.hpp"
#include "../random/rng.hpp"
#include "cvector.hpp"
#include <iostream>
#include <ranges>
namespace cp_algo::math::fft {
    template<modint_type base>
    struct dft {
        cvector A, B;
        static base factor, ifactor;
        static bool init;
        static int split;

        dft(auto const& a, size_t n): A(n), B(n) {
            if(!init) {
                factor = 1 + random::rng() % (base::mod() - 1);
                split = int(std::sqrt(base::mod())) + 1;
                ifactor = base(1) / factor;
                init = true;
            }
            base cur = factor;
            base step = bpow(factor, n);
            cvector::exec_on_roots(2 * n, std::min(n, size(a)), [&](size_t i, auto rt) {
                auto splt = [&](size_t i, auto mul) {
                    auto ai = i < size(a) ? (a[i] * mul).rem_direct() : 0;
                    auto rem = ai % split;
                    auto quo = (ai - rem) / split;
                    return std::pair{rem, quo};
                };
                auto [rai, qai] = splt(i, cur);
                auto [rani, qani] = splt(n + i, cur * step);
                A.set(i, point(rai, rani) * rt);
                B.set(i, point(qai, qani) * rt);
                cur *= factor;
            });
            checkpoint("dft init");
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
            for(size_t k = 0; k < n; k += flen) {
                auto rt = cvector::eval_point(k / flen / 2);
                if(k / flen % 2) {
                    rt = -rt;
                }
                auto [Ax, Ay] = A.at(k);
                auto [Bx, By] = B.at(k);
                vpoint AC, AD, BC, BD;
                AC = AD = BC = BD = vz;
                auto Cv = C.at(k), Dv = D.at(k);
                for (size_t i = 0; i < flen; i++) {
                    vpoint Av = {vz + Ax[i], vz + Ay[i]}, Bv = {vz + Bx[i], vz + By[i]};
                    AC += Av * Cv; AD += Av * Dv;
                    BC += Bv * Cv; BD += Bv * Dv;
                    real(Cv) = __builtin_shufflevector(real(Cv), real(Cv), 3, 0, 1, 2);
                    imag(Cv) = __builtin_shufflevector(imag(Cv), imag(Cv), 3, 0, 1, 2);
                    real(Dv) = __builtin_shufflevector(real(Dv), real(Dv), 3, 0, 1, 2);
                    imag(Dv) = __builtin_shufflevector(imag(Dv), imag(Dv), 3, 0, 1, 2);
                    auto cx = real(Cv)[0], cy = imag(Cv)[0];
                    auto dx = real(Dv)[0], dy = imag(Dv)[0];
                    real(Cv)[0] = cx * real(rt) - cy * imag(rt);
                    imag(Cv)[0] = cx * imag(rt) + cy * real(rt);
                    real(Dv)[0] = dx * real(rt) - dy * imag(rt);
                    imag(Dv)[0] = dx * imag(rt) + dy * real(rt);
                }
                A.at(k) = AC;
                C.at(k) = AD + BC;
                B.at(k) = BD;
            }
            checkpoint("dot");
            A.ifft();
            B.ifft();
            C.ifft();
            auto splitsplit = (base(split) * split).rem();
            base cur = ifactor * ifactor;
            base step = bpow(ifactor, n);
            cvector::exec_on_roots(2 * n, std::min(n, k), [&](size_t i, point rt) {
                rt = conj(rt);
                auto Ai = A.get(i) * rt;
                auto Bi = B.get(i) * rt;
                auto Ci = C.get(i) * rt;
                int64_t A0 = llround(real(Ai));
                int64_t A1 = llround(real(Ci));
                int64_t A2 = llround(real(Bi));
                res[i].setr_direct(base::m_reduce(A0 + A1 * split + A2 * splitsplit));
                res[i] *= cur;
                if(n + i < k) {
                    int64_t B0 = llround(imag(Ai));
                    int64_t B1 = llround(imag(Ci));
                    int64_t B2 = llround(imag(Bi));
                    res[n + i].setr_direct(base::m_reduce(B0 + B1 * split + B2 * splitsplit));
                    res[n + i] *= cur * step;
                }
                cur *= ifactor;
            });
            checkpoint("recover mod");
        }
        void mul_inplace(auto &&B, auto& res, size_t k) {
            mul(B.A, B.B, res, k);
        }
        void mul(auto const& B, auto& res, size_t k) {
            mul(cvector(B.A), B.B, res, k);
        }
        std::vector<base, big_alloc<base>> operator *= (dft &B) {
            std::vector<base, big_alloc<base>> res(2 * A.size());
            mul_inplace(B, res, size(res));
            return res;
        }
        std::vector<base, big_alloc<base>> operator *= (dft const& B) {
            std::vector<base, big_alloc<base>> res(2 * A.size());
            mul(B, res, size(res));
            return res;
        }
        auto operator * (dft const& B) const {
            return dft(*this) *= B;
        }

        point operator [](int i) const {return A.get(i);}
    };
    template<modint_type base> base dft<base>::factor = 1;
    template<modint_type base> base dft<base>::ifactor = 1;
    template<modint_type base> bool dft<base>::init = false;
    template<modint_type base> int dft<base>::split = 1;
    
    void mul_slow(auto &a, auto const& b, size_t k) {
        if(empty(a) || empty(b)) {
            a.clear();
        } else {
            size_t n = std::min(k, size(a));
            size_t m = std::min(k, size(b));
            a.resize(k);
            for(int j = int(k - 1); j >= 0; j--) {
                a[j] *= b[0];
                for(int i = std::max(j - (int)n, 0) + 1; i < std::min(j + 1, (int)m); i++) {
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
        if(std::min({k, size(a), size(b)}) < magic) {
            mul_slow(a, b, k);
            return;
        }
        auto n = std::max(flen, std::bit_ceil(
            std::min(k, size(a)) + std::min(k, size(b)) - 1
        ) / 2);
        auto A = dft<base>(a | std::views::take(k), n);
        a.assign(k, 0);
        checkpoint("reset a");
        if(&a == &b) {
            A.mul(A, a, k);
        } else {
            A.mul_inplace(dft<base>(b | std::views::take(k), n), a, k);
        }
    }
    void mul(auto &a, auto const& b) {
        size_t N = size(a) + size(b) - 1;
        if(std::max(size(a), size(b)) > (1 << 23)) {
            // do karatsuba to save memory
            auto n = (std::max(size(a), size(b)) + 1) / 2;
            auto a0 = to<std::vector>(a | std::views::take(n));
            auto a1 = to<std::vector>(a | std::views::drop(n));
            auto b0 = to<std::vector>(b | std::views::take(n));
            auto b1 = to<std::vector>(b | std::views::drop(n));
            a0.resize(n); a1.resize(n);
            b0.resize(n); b1.resize(n);
            auto a01 = to<std::vector>(std::views::zip_transform(std::plus{}, a0, a1));
            auto b01 = to<std::vector>(std::views::zip_transform(std::plus{}, b0, b1));
            mul(a0, b0);
            mul(a1, b1);
            mul(a01, b01);
            a.assign(4 * n, 0);
            for(auto [i, ai]: a0 | std::views::enumerate) {
                a[i] += ai;
                a[i + n] -= ai;
            }
            for(auto [i, ai]: a1 | std::views::enumerate) {
                a[i + n] -= ai;
                a[i + 2 * n] += ai;
            }
            for(auto [i, ai]: a01 | std::views::enumerate) {
                a[i + n] += ai;
            }
            a.resize(N);
        } else if(size(a)) {
            mul_truncate(a, b, N);
        }
    }
}
#endif // CP_ALGO_MATH_FFT_HPP
