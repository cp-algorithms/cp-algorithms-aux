#ifndef CP_ALGO_MATH_MULTIVAR_HPP
#define CP_ALGO_MATH_MULTIVAR_HPP
#include "../util/big_alloc.hpp"
#include "../number_theory/modint.hpp"
#include "../math/fft.hpp"
#pragma GCC push_options
#pragma GCC target("avx2")
namespace cp_algo::math::fft {
    template<modint_type base>
    struct multivar {
        big_vector<base> data;
        big_vector<size_t> ranks;
        big_vector<size_t> dim;
        size_t N;
        size_t rank(size_t i) {
            size_t cur = 1, res = 0, K = size(dim);
            for(auto ni: dim) {
                cur *= ni;
                res += i / cur;
            }
            return res % K;
        }
        multivar(auto const& dim): dim(begin(dim), end(dim)), N(
            std::ranges::fold_left(dim, 1, std::multiplies{})
        ) {
            data.resize(N);
            ranks.resize(N);
            for(auto [i, x]: ranks | std::views::enumerate) {
                x = rank(i);
            }
            checkpoint("multivar init");
        }
        void read() {
            for(auto &it: data) {
                std::cin >> it;
            }
            checkpoint("multivar read");
        }
        void print() {
            for(auto &it: data) {
                std::cout << it << " ";
            }
            std::cout << "\n";
            checkpoint("multivar write");
        }
        void mul(multivar<base> const& b) {
            assert(dim == b.dim);
            size_t K = size(dim);
            if(K == 0) {
                data[0] *= b.data[0];
                return;
            }
            big_vector<dft<base>> A, B;
            size_t M = std::max(flen, std::bit_ceil(2 * N - 1) / 2);
            for(size_t i = 0; i < K; i++) {
                A.emplace_back(data | std::views::enumerate | std::views::transform(
                    [&](auto jx) {
                        auto [j, x] = jx;
                        return ranks[j] == i ? x : base(0);
                    }
                ), M, false);
                B.emplace_back(b.data | std::views::enumerate | std::views::transform(
                    [&](auto jx) {
                        auto [j, x] = jx;
                        return ranks[j] == i ? x : base(0);
                    }
                ), M, false);
            }
            for(size_t i = 0; i < K; i++) {
                dft<base> C(M);
                cvector X = C.A;
                for(size_t j = 0; j < K; j++) {
                    size_t tj = (i - j + K) % K;
                    A[j].template dot<false, false>(B[tj].A, B[tj].B, C.A, C.B, X);
                }
                checkpoint("dot");
                big_vector<base> res((N + flen - 1) / flen * flen);
                C.A.template ifft<false>();
                C.B.template ifft<false>();
                X.template ifft<false>();
                C.recover_mod(X, res, N);
                for(size_t j = 0; j < N; j++) {
                    if(i == ranks[j]) {
                        data[j] = res[j];
                    }
                }
                checkpoint("store");
            }
        }
    };
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_MULTIVAR_HPP