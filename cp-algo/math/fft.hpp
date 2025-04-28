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
        using Int2 = base::Int2;
        static bool _init;
        static int split() {
            static const int splt = int(std::sqrt(base::mod())) + 1;
            return splt;
        }
        static u64x4 mod, imod;

        void init() {
            if(!_init) {
                factor = 1 + random::rng() % (base::mod() - 1);
                ifactor = base(1) / factor;
                mod = u64x4() + base::mod();
                imod = u64x4() + inv2(-base::mod());
                _init = true;
            }
        }

        dft(auto const& a, size_t n): A(n), B(n) {
            init();
            base b2x32 = bpow(base(2), 32);
            u64x4 cur = {
                (bpow(factor, 1) * b2x32).getr(),
                (bpow(factor, 2) * b2x32).getr(),
                (bpow(factor, 3) * b2x32).getr(),
                (bpow(factor, 4) * b2x32).getr()
            };
            u64x4 step4 = u64x4{} + (bpow(factor, 4) * b2x32).getr();
            u64x4 stepn = u64x4{} + (bpow(factor, n) * b2x32).getr();
            for(size_t i = 0; i < std::min(n, size(a)); i += flen) {
                auto splt = [&](size_t i, auto mul) {
                    if(i >= size(a)) {
                        return std::pair{vftype(), vftype()};
                    }
                    u64x4 au = {
                        i < size(a) ? a[i].getr() : 0,
                        i + 1 < size(a) ? a[i + 1].getr() : 0,
                        i + 2 < size(a) ? a[i + 2].getr() : 0,
                        i + 3 < size(a) ? a[i + 3].getr() : 0
                    };
                    au = montgomery_mul(au, mul, mod, imod);
                    au = au >= base::mod() ? au - base::mod() : au;
                    auto ai = to_double(i64x4(au >= base::mod() / 2 ? au - base::mod() : au));
                    auto quo = round(ai / split());
                    return std::pair{ai - quo * split(), quo};
                };
                auto [rai, qai] = splt(i, cur);
                auto [rani, qani] = splt(n + i, montgomery_mul(cur, stepn, mod, imod));
                A.at(i) = vpoint(rai, rani);
                B.at(i) = vpoint(qai, qani);
                cur = montgomery_mul(cur, step4, mod, imod);
            }
            checkpoint("dft init");
            if(n) {
                A.fft();
                B.fft();
            }
        }

        void dot(auto &&C, auto const& D) {
            cvector::exec_on_evals<1>(A.size() / flen, [&](size_t k, point rt) {
                k *= flen;
                auto [Ax, Ay] = A.at(k);
                auto [Bx, By] = B.at(k);
                vpoint AC, AD, BC, BD;
                AC = AD = BC = BD = vz;
                auto Cv = C.at(k), Dv = D.at(k);
                for (size_t i = 0; i < flen; i++) {
                    vpoint Av = {vz + Ax[i], vz + Ay[i]}, Bv = {vz + Bx[i], vz + By[i]};
                    AC += Av * Cv; AD += Av * Dv;
                    BC += Bv * Cv; BD += Bv * Dv;
                    real(Cv) = rotate_right(real(Cv));
                    imag(Cv) = rotate_right(imag(Cv));
                    real(Dv) = rotate_right(real(Dv));
                    imag(Dv) = rotate_right(imag(Dv));
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
            });
            checkpoint("dot");
        }

        void recover_mod(auto &&C, auto &res, size_t k) {
            res.assign((k / flen + 1) * flen, base(0));
            size_t n = A.size();
            auto const splitsplit = base(split() * split()).getr();
            base b2x32 = bpow(base(2), 32);
            base b2x64 = bpow(base(2), 64);
            u64x4 cur = {
                (bpow(ifactor, 2) * b2x64).getr(),
                (bpow(ifactor, 3) * b2x64).getr(),
                (bpow(ifactor, 4) * b2x64).getr(),
                (bpow(ifactor, 5) * b2x64).getr()
            };
            u64x4 step4 = u64x4{} + (bpow(ifactor, 4) * b2x32).getr();
            u64x4 stepn = u64x4{} + (bpow(ifactor, n) * b2x32).getr();
            for(size_t i = 0; i < std::min(n, k); i += flen) {
                auto [Ax, Ay] = A.at(i);
                auto [Bx, By] = B.at(i);
                auto [Cx, Cy] = C.at(i);
                auto set_i = [&](size_t i, auto A, auto B, auto C, auto mul) {
                    auto A0 = lround(A), A1 = lround(C), A2 = lround(B);
                    auto Ai = A0 + A1 * split() + A2 * splitsplit + uint64_t(base::modmod());
                    auto Au = montgomery_reduce(u64x4(Ai), mod, imod);
                    Au = montgomery_mul(Au, mul, mod, imod);
                    Au = Au >= base::mod() ? Au - base::mod() : Au;
                    for(size_t j = 0; j < flen; j++) {
                        res[i + j].setr(typename base::UInt(Au[j]));
                    }
                };
                set_i(i, Ax, Bx, Cx, cur);
                if(i + n < k) {
                    set_i(i + n, Ay, By, Cy, montgomery_mul(cur, stepn, mod, imod));
                }
                cur = montgomery_mul(cur, step4, mod, imod);
            }
            res.resize(k);
            checkpoint("recover mod");
        }

        void mul(auto &&C, auto const& D, auto &res, size_t k) {
            assert(A.size() == C.size());
            size_t n = A.size();
            if(!n) {
                res = {};
                return;
            }
            dot(C, D);
            A.ifft();
            B.ifft();
            C.ifft();
            recover_mod(C, res, k);
        }
        void mul_inplace(auto &&B, auto& res, size_t k) {
            mul(B.A, B.B, res, k);
        }
        void mul(auto const& B, auto& res, size_t k) {
            mul(cvector(B.A), B.B, res, k);
        }
        std::vector<base, big_alloc<base>> operator *= (dft &B) {
            std::vector<base, big_alloc<base>> res;
            mul_inplace(B, res, 2 * A.size());
            return res;
        }
        std::vector<base, big_alloc<base>> operator *= (dft const& B) {
            std::vector<base, big_alloc<base>> res;
            mul(B, res, 2 * A.size());
            return res;
        }
        auto operator * (dft const& B) const {
            return dft(*this) *= B;
        }

        point operator [](int i) const {return A.get(i);}
    };
    template<modint_type base> base dft<base>::factor = 1;
    template<modint_type base> base dft<base>::ifactor = 1;
    template<modint_type base> bool dft<base>::_init = false;
    template<modint_type base> u64x4 dft<base>::mod = {};
    template<modint_type base> u64x4 dft<base>::imod = {};
    
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
        if(&a == &b) {
            A.mul(A, a, k);
        } else {
            A.mul_inplace(dft<base>(b | std::views::take(k), n), a, k);
        }
    }
    void mul(auto &a, auto const& b) {
        size_t N = size(a) + size(b) - 1;
        if(std::max(size(a), size(b)) > (1 << 23)) {
            using T = std::decay_t<decltype(a[0])>;
            // do karatsuba to save memory
            auto n = (std::max(size(a), size(b)) + 1) / 2;
            auto a0 = to<std::vector<T, big_alloc<T>>>(a | std::views::take(n));
            auto a1 = to<std::vector<T, big_alloc<T>>>(a | std::views::drop(n));
            auto b0 = to<std::vector<T, big_alloc<T>>>(b | std::views::take(n));
            auto b1 = to<std::vector<T, big_alloc<T>>>(b | std::views::drop(n));
            a0.resize(n); a1.resize(n);
            b0.resize(n); b1.resize(n);
            auto a01 = to<std::vector<T, big_alloc<T>>>(std::views::zip_transform(std::plus{}, a0, a1));
            auto b01 = to<std::vector<T, big_alloc<T>>>(std::views::zip_transform(std::plus{}, b0, b1));
            checkpoint("karatsuba split");
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
            checkpoint("karatsuba join");
        } else if(size(a)) {
            mul_truncate(a, b, N);
        }
    }
}
#endif // CP_ALGO_MATH_FFT_HPP
