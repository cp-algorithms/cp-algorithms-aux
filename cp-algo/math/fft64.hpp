#ifndef CP_ALGO_MATH_FFT64_HPP
#define CP_ALGO_MATH_FFT64_HPP
#include "../random/rng.hpp"
#include "../math/common.hpp"
#include "../math/cvector.hpp"
#pragma GCC push_options
#pragma GCC target("avx2")

namespace cp_algo::math::fft {
    struct dft64 {
        big_vector<cp_algo::math::fft::cvector> cv;

        static uint64_t factor, ifactor;
        static bool _init;

        static void init() {
            if(_init) return;
            _init = true;
            factor = random::rng();
            if(factor % 2 == 0) {factor++;}
            ifactor = inv2(factor);
        }

        dft64(auto const& a, size_t n): cv(4, n) {
            init();
            uint64_t cur = 1, step = bpow(factor, n);
            for(size_t i = 0; i < std::min(std::size(a), n); i++) {
                auto split = [&](size_t i, uint64_t mul) -> std::array<int16_t, 4> {
                    uint64_t x = i < std::size(a) ? a[i] * mul : 0;
                    std::array<int16_t, 4> res;
                    for(int z = 0; z < 4; z++) {
                        res[z] = int16_t(x);
                        x = (x >> 16) + (res[z] < 0);
                    }
                    return res;
                };
                auto re = split(i, cur);
                auto im = split(n + i, cur * step);
                for(int z = 0; z < 4; z++) {
                    real(cv[z].at(i))[i % 4] = re[z];
                    imag(cv[z].at(i))[i % 4] = im[z];
                }
                cur *= factor;
            }
            checkpoint("dft64 init");
            for(auto &x: cv) {
                x.fft();
            }
        }

        static void do_dot_iter(point rt, std::array<vpoint, 4>& B, std::array<vpoint, 4> const& A, std::array<vpoint, 4>& C) {
            for(size_t k = 0; k < 4; k++) {
                for(size_t i = 0; i <= k; i++) {
                    C[k] += A[i] * B[k - i];
                }
            }
            for(size_t k = 0; k < 4; k++) {
                real(B[k]) = rotate_right(real(B[k]));
                imag(B[k]) = rotate_right(imag(B[k]));
                auto bx = real(B[k])[0], by = imag(B[k])[0];
                real(B[k])[0] = bx * real(rt) - by * imag(rt);
                imag(B[k])[0] = bx * imag(rt) + by * real(rt);
            }
        }

        void dot(dft64 const& t) {
            size_t N = cv[0].size();
            cvector::exec_on_evals<1>(N / flen, [&](size_t k, point rt) __attribute__((always_inline)) {
                k *= flen;
                auto [A0x, A0y] = cv[0].at(k);
                auto [A1x, A1y] = cv[1].at(k);
                auto [A2x, A2y] = cv[2].at(k);
                auto [A3x, A3y] = cv[3].at(k);
                std::array B = {
                    t.cv[0].at(k),
                    t.cv[1].at(k),
                    t.cv[2].at(k),
                    t.cv[3].at(k)
                };

                std::array<vpoint, 4> C = {vz, vz, vz, vz};
                for (size_t i = 0; i < flen; i++) {
                    std::array A = {
                        vpoint{vz + A0x[i], vz + A0y[i]},
                        vpoint{vz + A1x[i], vz + A1y[i]},
                        vpoint{vz + A2x[i], vz + A2y[i]},
                        vpoint{vz + A3x[i], vz + A3y[i]}
                    };
                    do_dot_iter(rt, B, A, C);
                }
                cv[0].at(k) = C[0];
                cv[1].at(k) = C[1];
                cv[2].at(k) = C[2];
                cv[3].at(k) = C[3];
            });
            checkpoint("dot");
            for(auto &x: cv) {
                x.ifft();
            }
        }

        void recover_mod(auto &res, size_t k) {
            size_t n = cv[0].size();
            uint64_t cur = 1, step = bpow(ifactor, n);
            for(size_t i = 0; i < std::min(k, n); i++) {
                std::array re = {real(cv[0].get(i)), real(cv[1].get(i)), real(cv[2].get(i)), real(cv[3].get(i))};
                std::array im = {imag(cv[0].get(i)), imag(cv[1].get(i)), imag(cv[2].get(i)), imag(cv[3].get(i))};
                auto set_i = [&](size_t i, auto &x, auto mul) {
                    if (i >= k) return;
                    res[i] = llround(x[0]) + (llround(x[1]) << 16) + (llround(x[2]) << 32) + (llround(x[3]) << 48);
                    res[i] *= mul;
                };
                set_i(i, re, cur);
                set_i(n + i, im, cur * step);
                cur *= ifactor;
            }
            cp_algo::checkpoint("recover mod");
        }
    };
    uint64_t dft64::factor = 1, dft64::ifactor = 1;
    bool dft64::_init = false;

    void conv64(auto& a, auto const& b) {
        size_t n = a.size(), m = b.size();
        size_t N = std::max(flen, std::bit_ceil(n + m - 1) / 2);
        dft64 A(a, N), B(b, N);
        A.dot(B);
        a.resize(n + m - 1);
        A.recover_mod(a, n + m - 1);
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_FFT64_HPP
