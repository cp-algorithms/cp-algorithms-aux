#ifndef CP_ALGO_ALGEBRA_POLY_IMPL_EUCLID_HPP
#define CP_ALGO_ALGEBRA_POLY_IMPL_EUCLID_HPP
#include "../../affine.hpp"
#include "../../fft.hpp"
#include <algorithm>
#include <numeric>
#include <cassert>
#include <vector>
#include <tuple>
#include <list>
// operations related to gcd and Euclidean algo
namespace cp_algo::algebra::poly::impl {
    template<typename poly>
    using gcd_result = std::pair<
        std::list<std::decay_t<poly>>,
        linfrac<std::decay_t<poly>>>;

    template<typename poly>
    gcd_result<poly> half_gcd(poly &&A, poly &&B) {
        assert(A.deg() >= B.deg());
        int m = size(A.a) / 2;
        if(B.deg() < m) {
            return {};
        }
        auto [ai, R] = A.divmod(B);
        std::tie(A, B) = {B, R};
        std::list a = {ai};
        auto T = -linfrac(ai).adj();

        auto advance = [&](int k) {
            auto [ak, Tk] = half_gcd(A.div_xk(k), B.div_xk(k));
            a.splice(end(a), ak);
            T.prepend(Tk);
            return Tk;
        };
        advance(m).apply(A, B);
        if constexpr (std::is_reference_v<poly>) {
            advance(2 * m - A.deg()).apply(A, B);
        } else {
            advance(2 * m - A.deg());
        }
        return {std::move(a), std::move(T)};
    }
    template<typename poly>
    gcd_result<poly> full_gcd(poly &&A, poly &&B) {
        using poly_t = std::decay_t<poly>;
        std::list<poly_t> ak;
        std::vector<linfrac<poly_t>> trs;
        while(!B.is_zero()) {
            auto [a0, R] = A.divmod(B);
            ak.push_back(a0);
            trs.push_back(-linfrac(a0).adj());
            std::tie(A, B) = {B, R};

            auto [a, Tr] = half_gcd(A, B);
            ak.splice(end(ak), a);
            trs.push_back(Tr);
        }
        return {ak, std::accumulate(rbegin(trs), rend(trs), linfrac<poly_t>{}, std::multiplies{})};
    }
}
#endif // CP_ALGO_ALGEBRA_POLY_IMPL_EUCLID_HPP
        