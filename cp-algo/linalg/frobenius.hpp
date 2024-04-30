#ifndef CP_ALGO_LINALG_FROBENIUS_HPP
#define CP_ALGO_LINALG_FROBENIUS_HPP
#include "matrix.hpp"
#include "../algebra/poly.hpp"
#include <algorithm>
#include <vector>
namespace cp_algo::linalg {
    enum frobenius_mode {blocks, full};
    template<frobenius_mode mode = blocks>
    auto frobenius_form(auto const& A) {
        using matrix = std::decay_t<decltype(A)>;
        using base = matrix::base;
        using polyn = algebra::poly_t<base>;
        assert(A.n() == A.m());
        size_t n = A.n();
        std::vector<polyn> charps;
        std::vector<vec<base>> basis, basis_init;
        while(size(basis) < n) {
            size_t start = size(basis);
            auto generate_block = [&](auto x) {
                while(true) {
                    vec<base> y = x | vec<base>::ei(n + 1, size(basis));
                    for(auto &it: basis) {
                        y.reduce_by(it);
                    }
                    y.normalize();
                    if(vec<base>(y[std::slice(0, n, 1)]) == vec<base>(n)) {
                        return polyn(std::vector<base>(begin(y) + n, end(y)));
                    } else {
                        basis_init.push_back(x);
                        basis.push_back(y);
                        x = A.apply(x);
                    }
                }
            };
            auto full_rec = generate_block(vec<base>::random(n));
            // Extra trimming to make it block-diagonal (expensive)
            if constexpr (mode == full) {
                if(full_rec.mod_xk(start) != polyn()) {
                    auto charp = full_rec.div_xk(start);
                    auto x = basis_init[start];
                    start = 0;
                    for(auto &rec: charps) {
                        polyn cur_rec = full_rec.substr(start, rec.deg());
                        auto shift = cur_rec / charp;
                        for(int j = 0; j <= shift.deg(); j++) {
                            x.add_scaled(basis_init[start + j], shift[j]);
                        }
                        start += rec.deg();
                    }
                    basis.resize(start);
                    basis_init.resize(start);
                    full_rec = generate_block(x.normalize());
                }
            }
            charps.push_back(full_rec.div_xk(start));
        }
        // Find transform matrices while we're at it...
        if constexpr (mode == full) {
            for(size_t i = 0; i < size(basis); i++) {
                for(size_t j = i + 1; j < size(basis); j++) {
                    basis[i].reduce_by(basis[j]);
                }
                basis[i].normalize();
                basis[i] = vec<base>(
                    basis[i][std::slice(n, n, 1)]
                ) * (base(1) / basis[i][i]);
            }
            auto T = matrix::from_range(basis_init);
            auto Tinv = matrix::from_range(basis);
            return std::tuple{T, Tinv, charps};
        } else {
            return charps;
        }
    }

    template<typename base>
    auto frobenius_pow(matrix<base> A, uint64_t k) {
        using polyn = algebra::poly_t<base>;
        auto [T, Tinv, charps] = frobenius_form<full>(A);
        std::vector<matrix<base>> blocks;
        for(auto charp: charps) {
            matrix<base> block(charp.deg());
            auto xk = polyn::xk(1).powmod(k, charp);
            for(size_t i = 0; i < block.n(); i++) {
                std::ranges::copy(xk.a, begin(block[i]));
                xk = xk.mul_xk(1) % charp;
            }
            blocks.push_back(block);
        }
        auto S = matrix<base>::block_diagonal(blocks);
        return Tinv * S * T;
    }
};
#endif // CP_ALGO_LINALG_FROBENIUS_HPP
