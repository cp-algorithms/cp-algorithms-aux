#ifndef CP_ALGO_MATH_FFT_HPP
#define CP_ALGO_MATH_FFT_HPP
#include "../number_theory/modint.hpp"
#include "cvector.hpp"
namespace cp_algo::math::fft {
    template<typename base>
    struct dft {
        ftvec A;
        
        dft(std::vector<base> const& a, size_t n): A(n) {
            for(size_t i = 0; i < std::min(n, a.size()); i++) {
                A[i] = a[i];
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
                res[k] = A[k];
            }
            return res;
        }

        auto operator * (dft const& B) const {
            return dft(*this) *= B;
        }
    };

    template<modint_type base>
    struct dft<base> {
        int split;
        ftvec A, B;
        
        dft(auto const& a, size_t n): A(n), B(n) {
            n = size(A);
            split = int(std::sqrt(base::mod())) + 1;
            ftvec::exec_on_roots(2 * n, size(a), [&](size_t i, point rt) {
                size_t ti = std::min(i, i - n);
                auto rem = std::remainder(a[i].rem(), split);
                auto quo = (a[i].rem() - rem) / split;
                A[ti] += rem * rt;
                B[ti] += quo * rt;
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
            for(size_t i = 0; i < n; i += ftvec::threshold) {
                auto AC = ftvec::dot_block(i, A, C);
                auto AD = ftvec::dot_block(i, A, D);
                auto BC = ftvec::dot_block(i, B, C);
                auto BD = ftvec::dot_block(i, B, D);
                for(size_t j = 0; j < ftvec::threshold; j++) {
                    A[i + j] = AC[j];
                    C[i + j] = AD[j] + BC[j];
                    B[i + j] = BD[j];
                }
            }
            A.ifft();
            B.ifft();
            C.ifft();
            auto splitsplit = (base(split) * split).rem();
            ftvec::exec_on_roots(2 * n, std::min(n, k), [&](size_t i, point rt) {
                rt = conj(rt);
                auto Ai = A[i] * rt;
                auto Bi = B[i] * rt;
                auto Ci = C[i] * rt;
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
            mul(ftvec(B.A), B.B, res, k);
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
        return std::bit_ceil(as + bs - 1) / 2;
    }
    void mul_truncate(auto &a, auto const& b, size_t k) {
        using base = std::decay_t<decltype(a[0])>;
        if(std::min({k, size(a), size(b)}) < 1) {
            mul_slow(a, b, k);
            return;
        }
        auto n = std::bit_ceil(std::min(k, size(a)) + std::min(k, size(b)) - 1) / 2;
        a.resize(k);
        auto A = dft<base>(a, n);
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
