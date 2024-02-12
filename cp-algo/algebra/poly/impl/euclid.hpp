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
    auto half_gcd(poly &&A, poly &&B) {
        assert(A.deg() >= B.deg());
        using poly_t = std::decay_t<poly>;
        std::list<poly_t> a;
        linfrac<poly_t> T;
        int m = size(A.a) / 2;
        if(B.deg() < m) {
            return std::pair{a, T};
        }
        auto [ar, Tr] = half_gcd(A.div_xk(m), B.div_xk(m));
        std::tie(A, B) = Tr.apply(A, B);
        a = ar;
        T = Tr;
        if(B.deg() < m) {
            return std::pair{a, T};
        }
        auto [ai, R] = A.divmod(B);
        std::tie(A, B) = {B, R};
        int k = 2 * m - A.deg();
        auto [as, Ts] = half_gcd(A.div_xk(k), B.div_xk(k));
        a.push_back(ai);
        a.splice(end(a), as);
        T = Ts * -linfrac(ai).adj() * T;
        if constexpr (std::is_reference_v<poly>) {
            std::tie(A, B) = Ts.apply(A, B);
        }
        return std::pair{a, T};
    }
    template<typename poly>
    auto full_gcd(poly &&A, poly &&B) {
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
        std::ranges::reverse(trs);
        return std::pair{ak, std::accumulate(begin(trs), end(trs), linfrac<poly_t>{}, std::multiplies{})};
    }
}
#endif // CP_ALGO_ALGEBRA_POLY_IMPL_EUCLID_HPP
        