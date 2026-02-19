#ifndef CP_ALGO_MATH_MULTIVAR_HPP
#define CP_ALGO_MATH_MULTIVAR_HPP
#include "../util/big_alloc.hpp"
#include "../number_theory/modint.hpp"
#include "../math/fft.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::fft {
    template<modint_type base>
    struct multivar {
        big_vector<base> data;
        big_vector<size_t> ranks;
        big_vector<size_t> dim;
        size_t N;
        size_t rank(size_t i) {
            size_t cur = 1, res = 0, K = size(dim);
            if(K == 0) return 0;
            for(auto ni: dim) {
                cur *= ni;
                res += i / cur;
            }
            return res % K;
        }
        static void copy_prefix(
            big_vector<base>& dst,
            big_vector<size_t> const& dst_dim,
            big_vector<base> const& src,
            big_vector<size_t> const& src_dim,
            big_vector<size_t> const& iter_dim,
            size_t iter_N
        ) {
            size_t K = iter_dim.size();
            if(K == 0) {
                dst[0] = src[0];
                return;
            }
            if(K == 1) {
                std::copy_n(src.data(), iter_dim[0], dst.data());
                return;
            }
            if(K == 2) {
                size_t rows = iter_dim[1], cols = iter_dim[0];
                for(size_t j = 0; j < rows; j++) {
                    std::copy_n(
                        src.data() + j * src_dim[0],
                        cols,
                        dst.data() + j * dst_dim[0]
                    );
                }
                return;
            }
            big_vector<size_t> src_stride(K), dst_stride(K);
            src_stride[0] = 1;
            dst_stride[0] = 1;
            for(size_t i = 1; i < K; i++) {
                src_stride[i] = src_stride[i - 1] * src_dim[i - 1];
                dst_stride[i] = dst_stride[i - 1] * dst_dim[i - 1];
            }
            big_vector<size_t> idx(K);
            size_t src_index = 0, dst_index = 0;
            for(size_t t = 0; t < iter_N; t++) {
                dst[dst_index] = src[src_index];
                for(size_t d = 0; d < K; d++) {
                    idx[d]++;
                    src_index += src_stride[d];
                    dst_index += dst_stride[d];
                    if(idx[d] < iter_dim[d]) {
                        break;
                    }
                    idx[d] = 0;
                    src_index -= src_stride[d] * iter_dim[d];
                    dst_index -= dst_stride[d] * iter_dim[d];
                }
            }
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
        size_t linear_index(auto const& idx) const {
            size_t pos = 0, stride = 1;
            size_t K = dim.size();
            for(size_t i = 0; i < K; i++) {
                pos += idx[i] * stride;
                stride *= dim[i];
            }
            return pos;
        }
        template<class Idx> requires requires(Idx const& idx) { idx[0]; }
        base& operator[](Idx const& idx) {
            return data[linear_index(idx)];
        }
        template<class Idx> requires requires(Idx const& idx) { idx[0]; }
        base const& operator[](Idx const& idx) const {
            return data[linear_index(idx)];
        }
        template<std::convertible_to<size_t>... Args>
        base& operator[](Args... args) {
            size_t idx[] = {static_cast<size_t>(args)...};
            return data[linear_index(idx)];
        }
        template<std::convertible_to<size_t>... Args>
        base const& operator[](Args... args) const {
            size_t idx[] = {static_cast<size_t>(args)...};
            return data[linear_index(idx)];
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
        void assign_prefix_from(multivar<base> const& src) {
            assert(dim.size() == src.dim.size());
            size_t K = dim.size();
            if(K == 0) {
                data[0] = src.data[0];
                return;
            }
            for(size_t i = 0; i < K; i++) {
                assert(src.dim[i] <= dim[i]);
            }
            copy_prefix(data, dim, src.data, src.dim, src.dim, src.N);
        }
        multivar<base> truncated(auto const& new_dim) const {
            big_vector<size_t> nd(begin(new_dim), end(new_dim));
            assert(nd.size() == dim.size());
            for(size_t i = 0; i < nd.size(); i++) {
                assert(nd[i] <= dim[i]);
            }
            multivar<base> out(nd);
            if(out.N == 0) {
                return out;
            }
            copy_prefix(out.data, out.dim, data, dim, out.dim, out.N);
            return out;
        }
        void truncate_inplace(auto const& new_dim) {
            big_vector<size_t> nd(begin(new_dim), end(new_dim));
            assert(nd.size() == dim.size());
            size_t K = nd.size();
            for(size_t i = 0; i < K; i++) {
                assert(nd[i] <= dim[i]);
            }
            size_t new_N = std::ranges::fold_left(nd, 1, std::multiplies{});
            if(new_N == 0) {
                data.clear();
                dim = std::move(nd);
                ranks.clear();
                N = 0;
                return;
            }
            // In-place: copy elements forward (src >= dst, so no overlap issues)
            if(K == 1) {
                // 1D: just resize
            } else if(K == 2) {
                size_t rows = nd[1], cols = nd[0];
                for(size_t j = 1; j < rows; j++) {
                    std::copy_n(
                        data.data() + j * dim[0],
                        cols,
                        data.data() + j * cols
                    );
                }
            } else {
                copy_prefix(data, nd, data, dim, nd, new_N);
            }
            data.resize(new_N);
            dim = std::move(nd);
            N = new_N;
            ranks.resize(N);
            for(auto [i, x]: ranks | std::views::enumerate) {
                x = rank(i);
            }
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