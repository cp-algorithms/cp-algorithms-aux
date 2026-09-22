#ifndef CP_ALGO_LINALG_BLOCK_INVERSE_HPP
#define CP_ALGO_LINALG_BLOCK_INVERSE_HPP
#include "matrix.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::linalg::impl {
    // Shoup multiplication: for canonical inputs, product - quotient * mod < 2 * mod.
    // Keep 32-bit rows canonical after every update, restricted to [first, last).
    template<class base>
    void inverse_add_scaled(modint_vec<base> &x, modint_vec<base> const& y,
                            base factor, size_t first, size_t last) {
        constexpr uint32_t mod = base::mod();
        uint32_t scale = factor.getr();
        if(!scale) return;
        uint32_t quotient = (uint64_t(scale) << 32) / mod;
        for(; first + 8 <= last; first += 8) {
            u32x8 value, old;
            std::memcpy(&value, y.data() + first, sizeof value);
            std::memcpy(&old, x.data() + first, sizeof old);
            auto q = u64x4(u32x8() + quotient);
            auto lo = strassen_product<mod>::mul(u64x4(value), q) >> 32;
            auto hi = (strassen_product<mod>::mul(u64x4(value) >> 32, q) >> 32) << 32;
            auto product = value * scale - u32x8(lo | hi) * mod;
            auto result = reduce_once(old + reduce_once(product, mod), mod);
            std::memcpy(x.data() + first, &result, sizeof result);
        }
        for(; first < last; first++) x[first] += factor * y[first];
    }
}
namespace cp_algo::linalg {
    // Schur-complement inversion, suited to large matrices with invertible leading blocks.
    // Singular leading blocks fall back to Gaussian elimination, which can be faster there.
    template<class base, class row>
    std::pair<base, matrix<base, row>> block_inverse(matrix<base, row> const& A) {
        assert(A.n() == A.m());
        using matrix = linalg::matrix<base, row>;
        if constexpr(!impl::use_strassen<row>) {
            return A.inv();
        } else if constexpr(base::bits > 32) {
            // Convert once around the entire recursion; narrow callers avoid both copies.
            using small = math::modint<int(base::mod())>;
            linalg::matrix<small> a(A.n());
            for(size_t i = 0; i < A.n(); i++) for(size_t j = 0; j < A.n(); j++) {
                a[i][j].setr(A[i][j].getr() % base::mod());
            }
            auto [det, inverse] = block_inverse(a);
            matrix result(inverse.n());
            for(size_t i = 0; i < inverse.n(); i++) for(size_t j = 0; j < inverse.n(); j++) {
                result[i][j].setr(inverse[i][j].getr());
            }
            return {base(det.getr()), std::move(result)};
        } else {
            if(A.n() >= 128) {
                auto lo = std::views::take(A.n() / 2);
                auto hi = std::views::drop(A.n() / 2);
                auto [da, ai] = block_inverse(matrix(A.submatrix(lo, lo)));
                if(da != base(0)) {
                    matrix b = A.submatrix(lo, hi), c = A.submatrix(hi, lo), d = A.submatrix(hi, hi);
                    b.normalize(); c.normalize(); d.normalize();
                    auto mul = [](matrix const& x, matrix const& y) {
                        return impl::strassen_product<base::mod()>::product(x, y);
                    };
                    auto v = mul(c, ai);
                    d -= mul(v, b); // Schur complement D - C A^-1 B.
                    auto [ds, si] = block_inverse(d);
                    if(ds == base(0)) return {0, {}};
                    auto u = mul(ai, b);
                    auto r = mul(u, si), t = mul(si, v);
                    ai += mul(r, v);
                    matrix res(A.n());
                    size_t cut = ai.n();
                    for(size_t i = 0; i < cut; i++) {
                        std::ranges::copy(ai[i], res[i].begin());
                        for(size_t j = 0; j < r.m(); j++) res[i][cut + j] = -r[i][j];
                    }
                    for(size_t i = 0; i < si.n(); i++) {
                        for(size_t j = 0; j < cut; j++) res[cut + i][j] = -t[i][j];
                        std::ranges::copy(si[i], res[cut + i].begin() + cut);
                    }
                    return {da * ds, std::move(res)};
                }
            }
            auto add_scaled = impl::inverse_add_scaled<base>;
            size_t n = A.n();
            matrix a(A.submatrix(std::views::all, std::views::all));
            base det = 1;
            big_vector<size_t> permutation(n);
            for(size_t i = 0; i < n; i++) permutation[i] = i;
            // Factor P*A = L*U; store reciprocal pivots and a unit-diagonal U.
            for(size_t i = 0; i < n; i++) {
                size_t p = i;
                while(p < n && a[p][i] == base(0)) p++;
                if(p == n) return {0, {}};
                if(p != i) {
                    std::swap(a[p], a[i]); std::swap(permutation[p], permutation[i]);
                    det = -det;
                }
                det *= a[i][i];
                base inv = a[i][i].inv();
                // x *= inv via x += (inv - 1) * x, using the vectorized update.
                add_scaled(a[i], a[i], inv - base(1), i + 1, n);
                a[i][i] = inv;
                for(size_t j = i + 1; j < n; j++) if(a[j][i] != base(0)) {
                    add_scaled(a[j], a[i], -a[j][i], i + 1, n);
                }
            }
            // Delay inverse work until all pivots succeed, then solve L*Y = P.
            matrix inverse(n);
            big_vector<std::pair<size_t, size_t>> ranges(n);
            for(size_t i = 0; i < n; i++) {
                inverse[i][permutation[i]] = 1;
                ranges[i] = {permutation[i], permutation[i] + 1};
            }
            for(size_t i = 0; i < n; i++) {
                auto [left, right] = ranges[i];
                add_scaled(inverse[i], inverse[i], a[i][i] - base(1), left, right);
                for(size_t j = i + 1; j < n; j++) if(a[j][i] != base(0)) {
                    add_scaled(inverse[j], inverse[i], -a[j][i], left, right);
                    ranges[j].first = std::min(ranges[j].first, left);
                    ranges[j].second = std::max(ranges[j].second, right);
                }
            }
            // Back substitution in U*X = Y; U has unit diagonal.
            for(size_t i = n; i-- > 0;) {
                auto [left, right] = ranges[i];
                for(size_t j = 0; j < i; j++) if(a[j][i] != base(0)) {
                    add_scaled(inverse[j], inverse[i], -a[j][i], left, right);
                    ranges[j].first = std::min(ranges[j].first, left);
                    ranges[j].second = std::max(ranges[j].second, right);
                }
            }
            return {det, std::move(inverse)};
        }
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_LINALG_BLOCK_INVERSE_HPP
