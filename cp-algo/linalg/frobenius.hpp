#ifndef CP_ALGO_LINALG_FROBENIUS_HPP
#define CP_ALGO_LINALG_FROBENIUS_HPP
#include "../math/poly/powmod.hpp"
#include "matrix.hpp"
#include <algorithm>
#include <vector>
#include <ranges>
namespace cp_algo::linalg {
    enum frobenius_mode {blocks, full};
    template<frobenius_mode mode = blocks>
    auto frobenius_form(auto const& A) {
        using matrix = std::decay_t<decltype(A)>;
        using vec_t = matrix::vec_t;
        using base = matrix::base;
        using polyn = math::poly_t<base>;
        assert(A.n() == A.m());
        size_t n = A.n();
        big_vector<polyn> charps;
        big_vector<vec_t> basis, basis_init;
        while(size(basis) < n) {
            size_t start = size(basis);
            auto generate_block = [&](auto x) {
                while(true) {
                    vec_t y(2 * n + 1);
                    std::ranges::copy(x, begin(y));
                    y[n + size(basis)] = 1;
                    for(auto &it: basis) {
                        y.reduce_by(it);
                    }
                    y.normalize();
                    if(std::ranges::count(y | std::views::take(n), base(0)) == int(n)) {
                        return polyn(typename polyn::Vector(begin(y) + n, end(y)));
                    } else {
                        basis.push_back(std::move(y));
                        if constexpr(mode == full) {
                            basis_init.push_back(std::move(x));
                            x = A.apply(basis_init.back());
                        } else {
                            x = A.apply(x);
                        }
                    }
                }
            };
            auto full_rec = generate_block(vec_t::random(n));
            // Extra trimming to make it block-diagonal (expensive)
            if constexpr (mode == full) {
                if(full_rec.mod_xk(start) != polyn()) {
                    auto charp = full_rec.div_xk(start);
                    auto x = basis_init[start];
                    auto shift = full_rec / charp;
                    for(int j = 0; j < shift.deg(); j++) {
                        x.add_scaled(basis_init[j], shift[j]);
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
            for(size_t i = 0; i < n; i++) {
                for(size_t j = i + 1; j < n; j++) {
                    basis[i].reduce_by(basis[j]);
                }
                basis[i].normalize();
            }
            auto T = matrix(std::move(basis_init));
            auto Tinv = matrix(std::move(basis));
            std::ignore = Tinv.sort_classify(n);
            for(size_t i = 0; i < n; i++) {
                Tinv[i] = vec_t(
                    Tinv[i] | std::views::drop(n) | std::views::take(n)
                ) * (base(1) / Tinv[i][i]);
            }
            return std::tuple{std::move(T), std::move(Tinv), std::move(charps)};
        } else {
            return charps;
        }
    }

    template<typename base>
    auto with_frobenius(matrix<base> const& A, auto &&callback) {
        auto [T, Tinv, charps] = frobenius_form<full>(A);
        // Apply each diagonal block to its rows of T without forming S.
        matrix<base> ST(A.n());
        size_t start = 0;
        for(auto const& charp: charps) {
            size_t d = charp.deg();
            auto xk = callback(charp);
            for(size_t i = 0; i < d; i++) {
                for(size_t j = 0; j < xk.a.size(); j++) {
                    ST[start + i].add_scaled(T[start + j], xk[j]);
                }
                if(i + 1 < d) {
                    xk = xk.mul_xk(1) % charp;
                }
            }
            start += d;
        }
        return Tinv * ST.normalize();
    }

    template<typename base>
    auto frobenius_pow(matrix<base> const& A, uint64_t k) {
        return with_frobenius(A, [k](auto const& charp) {
            if(charp.deg() == 1) {
                return math::poly_t<base>(bpow(-charp[0] / charp[1], k));
            }
            return math::powmod(math::poly_t<base>::xk(1), k, charp);
        });
    }
};
#endif // CP_ALGO_LINALG_FROBENIUS_HPP
