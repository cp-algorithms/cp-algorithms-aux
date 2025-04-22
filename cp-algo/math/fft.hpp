#ifndef CP_ALGO_MATH_FFT_HPP
#define CP_ALGO_MATH_FFT_HPP
#include "../number_theory/modint.hpp"
#include "../util/checkpoint.hpp"
#include "cvector.hpp"
#include <ranges>
#include <iostream>
namespace cp_algo::math::fft {
    template<modint_type base>
    struct dft {
        int split;
        cvector A, B;
        
        dft(auto const& a, size_t n): A(n), B(n) {
            split = int(std::sqrt(base::mod())) + 1;
            cvector::exec_on_roots(2 * n, std::min(n, size(a)), [&](size_t i, auto rt) {
                auto splt = [&](size_t i) {
                    auto ai = ftype(i < size(a) ? a[i].rem() : 0);
                    auto rem = std::remainder(ai, split);
                    auto quo = (ai - rem) / split;
                    return std::pair{rem, quo};
                };
                auto [rai, qai] = splt(i);
                auto [rani, qani] = splt(n + i);
                A.set(i, point(rai, rani) * rt);
                B.set(i, point(qai, qani) * rt);
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
                auto [Ax, Ay] = A.vget(k);
                auto [Bx, By] = B.vget(k);
                auto [Cvx, Cvy] = C.vget(k);
                auto [Dvx, Dvy] = D.vget(k);
                auto [Crvx, Crvy] = vpoint(Cvx, Cvy) * vpoint(fz + real(rt), fz + imag(rt));
                auto [Drvx, Drvy] = vpoint(Dvx, Dvy) * vpoint(fz + real(rt), fz + imag(rt));
                vftype Cx[2] = {Crvx, Cvx}, Cy[2] = {Crvy, Cvy};
                vftype Dx[2] = {Drvx, Dvx}, Dy[2] = {Drvy, Dvy};
                vpoint AC, AD, BC, BD;
                AC = AD = BC = BD = {fz, fz};
                for(size_t i = 0; i < flen; i++) {
                    auto Csx = (vftype*)((ftype*)Cx + flen - i);
                    auto Csy = (vftype*)((ftype*)Cy + flen - i);
                    auto Dsx = (vftype*)((ftype*)Dx + flen - i);
                    auto Dsy = (vftype*)((ftype*)Dy + flen - i);
                    vpoint As = {fz + Ax[i], fz + Ay[i]}, Bs = {fz + Bx[i], fz + By[i]};
                    vpoint Cs = {*Csx, *Csy}, Ds = {*Dsx, *Dsy};
                    AC += As * Cs; AD += As * Ds;
                    BC += Bs * Cs; BD += Bs * Ds;
                }
                A.set(k, AC);
                C.set(k, AD + BC);
                B.set(k, BD);
            }
            checkpoint("dot");
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
            checkpoint("recover mod");
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
        a.resize(k);
        checkpoint("resize a");
        if(&a == &b) {
            A.mul(A, a, k);
        } else {
            A.mul_inplace(dft<base>(b | std::views::take(k), n), a, k);
        }
    }
    void mul(auto &a, auto const& b) {
        if(size(a)) {
            mul_truncate(a, b, size(a) + size(b) - 1);
        }
    }
}
#endif // CP_ALGO_MATH_FFT_HPP
