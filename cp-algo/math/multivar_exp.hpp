#ifndef CP_ALGO_MATH_MULTIVAR_EXP_HPP
#define CP_ALGO_MATH_MULTIVAR_EXP_HPP
#include "multivar_inv.hpp"
#include <algorithm>
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::fft {
    // Compute logarithm of multivariate formal power series
    // Requires: a.data[0] == 1
    // Returns: log(a) mod truncation
    //
    // Uses: d/dx_0 log(a) = (d/dx_0 a) / a, then integrate,
    // with the constant of integration = log(a|_{x_0=0}).
    template<modint_type base>
    multivar<base> multivar_log(multivar<base> const& a, auto const& target_dim) {
        assert(a.data[0] == base(1));

        size_t K = a.dim.size();
        big_vector<size_t> dim_vec(K);
        for(size_t i = 0; i < K; i++) {
            dim_vec[i] = target_dim[i];
        }

        if(K == 0) {
            return multivar<base>(dim_vec);
        }

        auto a_trunc = a.truncated(dim_vec);

        multivar<base> a_deriv(dim_vec);
        for(size_t j = 0; j < a_trunc.N; j++) {
            size_t idx0 = j % a_trunc.dim[0];
            if(idx0 > 0) {
                a_deriv.data[j - 1] += a_trunc.data[j] * base(idx0);
            }
        }

        auto a_inv = multivar_inv(a_trunc, dim_vec);

        a_deriv.mul(a_inv);

        multivar<base> result(dim_vec);
        for(size_t j = 0; j < a_deriv.N; j++) {
            size_t idx0 = j % dim_vec[0];
            if(idx0 + 1 < dim_vec[0]) {
                result.data[j + 1] = a_deriv.data[j] / base(idx0 + 1);
            }
        }

        big_vector<size_t> slice_dim(K - 1);
        for(size_t i = 1; i < K; i++) {
            slice_dim[i - 1] = dim_vec[i];
        }
        size_t slice_N = std::ranges::fold_left(slice_dim, 1, std::multiplies{});
        multivar<base> a_slice(slice_dim);
        for(size_t j = 0; j < slice_N; j++) {
            a_slice.data[j] = a_trunc.data[j * dim_vec[0]];
        }

        auto log_slice = multivar_log(a_slice, slice_dim);

        for(size_t j = 0; j < slice_N; j++) {
            result.data[j * dim_vec[0]] += log_slice.data[j];
        }

        return result;
    }

    // Compute exponential of multivariate formal power series using Newton iteration
    // Requires: a.data[0] == 0
    // Returns: exp(a) mod truncation
    //
    // Carries b_inv across iterations to avoid recomputing inv from scratch
    // in each log call. Per Newton step: 5 muls instead of 6.
    template<modint_type base>
    multivar<base> multivar_exp(multivar<base> const& a, auto const& target_dim) {
        assert(a.data[0] == base(0));

        size_t K = a.dim.size();

        big_vector<size_t> init_dim(K, 1);
        multivar<base> b(init_dim);
        b.data[0] = base(1);
        multivar<base> b_inv(init_dim);
        b_inv.data[0] = base(1);

        size_t degree = 1;
        size_t total_degree = 0;
        for(size_t i = 0; i < K; i++) {
            total_degree += target_dim[i];
        }

        big_vector<size_t> next_dim(K);

        while(degree < total_degree) {
            degree *= 2;

            for(size_t i = 0; i < K; i++) {
                next_dim[i] = std::min(target_dim[i], degree);
            }

            {
                multivar<base> b_ext(next_dim);
                b_ext.assign_prefix_from(b);
                b = std::move(b_ext);
            }
            {
                multivar<base> bi_ext(next_dim);
                bi_ext.assign_prefix_from(b_inv);
                b_inv = std::move(bi_ext);
            }

            // Refine b_inv: b_inv = b_inv * (2 - b * b_inv)
            {
                auto t = b;
                t.mul(b_inv);
                for(auto &x: t.data) x = -x;
                t.data[0] += base(2);
                b_inv.mul(t);
            }

            // Compute log(b) using the carried b_inv
            multivar<base> log_b(next_dim);
            {
                multivar<base> b_deriv(next_dim);
                for(size_t j = 0; j < b.N; j++) {
                    size_t idx0 = j % b.dim[0];
                    if(idx0 > 0) {
                        b_deriv.data[j - 1] += b.data[j] * base(idx0);
                    }
                }

                b_deriv.mul(b_inv);

                for(size_t j = 0; j < b_deriv.N; j++) {
                    size_t idx0 = j % next_dim[0];
                    if(idx0 + 1 < next_dim[0]) {
                        log_b.data[j + 1] = b_deriv.data[j] / base(idx0 + 1);
                    }
                }

                if(K > 1) {
                    big_vector<size_t> slice_dim(K - 1);
                    for(size_t i = 1; i < K; i++) {
                        slice_dim[i - 1] = next_dim[i];
                    }
                    size_t slice_N = std::ranges::fold_left(
                        slice_dim, 1, std::multiplies{}
                    );
                    multivar<base> b_slice(slice_dim);
                    for(size_t j = 0; j < slice_N; j++) {
                        b_slice.data[j] = b.data[j * next_dim[0]];
                    }
                    auto log_slice = multivar_log(b_slice, slice_dim);
                    for(size_t j = 0; j < slice_N; j++) {
                        log_b.data[j * next_dim[0]] += log_slice.data[j];
                    }
                }
            }

            // eps = a - log(b)
            auto c = a.truncated(next_dim);
            for(size_t i = 0; i < c.N; i++) {
                c.data[i] -= log_b.data[i];
            }

            // b *= (1 + eps)
            c.data[0] += base(1);
            b.mul(c);

            // b_inv *= (1 - eps)
            for(auto &x: c.data) x = -x;
            c.data[0] += base(2);
            b_inv.mul(c);
        }

        b.truncate_inplace(target_dim);
        return b;
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_MULTIVAR_EXP_HPP
