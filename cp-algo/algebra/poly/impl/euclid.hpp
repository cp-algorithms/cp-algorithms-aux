#ifndef CP_ALGO_ALGEBRA_POLY_IMPL_EUCLID_HPP
#define CP_ALGO_ALGEBRA_POLY_IMPL_EUCLID_HPP
#include "../../affine.hpp"
#include <cassert>
#include <vector>
#include <tuple>
// operations related to gcd and Euclidean algo
namespace cp_algo::algebra::poly::impl {
    void concat(auto &a, auto const& b) {
        a.insert(a.end(), b.begin(), b.end());
    }
    template<typename poly>
    std::pair<std::vector<poly>, linfrac<poly>> half_gcd(poly A, poly B) {
        assert(A.deg() >= B.deg());
        int m = (A.deg() + 1) / 2;
        if(B.deg() < m) {
            return {};
        }
        auto [ar, Tr] = half_gcd(A.div_xk(m), B.div_xk(m));
        std::tie(A, B) = Tr.adj().apply(A, B);
        if(B.deg() < m) {
            return std::pair{ar, Tr};
        }
        auto [ai, R] = A.divmod(B);
        std::tie(A, B) = {B, R};
        int k = 2 * m - A.deg();
        auto [as, Ts] = half_gcd(A.div_xk(k), B.div_xk(k));
        ar.push_back(ai);
        concat(ar, as);
        return std::pair{ar, Tr * linfrac(ai) * Ts};
    }
    template<typename poly>
    auto full_gcd(poly A, poly B) {
        std::vector<poly> ak;
        std::vector<linfrac<poly>> trs;
        while(!B.is_zero()) {
            if(2 * B.deg() > A.deg()) {
                auto [a, Tr] = half_gcd(A, B);
                concat(ak, a);
                trs.push_back(Tr);
                std::tie(A, B) = trs.back().adj().apply(A, B);
            } else {
                auto [a, R] = A.divmod(B);
                ak.push_back(a);
                trs.emplace_back(a);
                std::tie(A, B) = {B, R};
            }
        }
        trs.emplace_back();
        while(trs.size() >= 2) {
            trs[trs.size() - 2] = trs[trs.size() - 2] * trs[trs.size() - 1];
            trs.pop_back();
        }
        return std::pair{ak, trs.back()};
    }
}
#endif // CP_ALGO_ALGEBRA_POLY_IMPL_EUCLID_HPP
        