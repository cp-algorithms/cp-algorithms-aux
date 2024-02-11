#ifndef CP_ALGO_ALGEBRA_POLY_IMPL_DIV_HPP
#define CP_ALGO_ALGEBRA_POLY_IMPL_DIV_HPP
#include "../../common.hpp"
#include <cassert>
// operations related to polynomial division
namespace cp_algo::algebra::poly::impl {
    auto divmod_slow(auto const& p, auto const& q) {
        auto R = p;
        decltype(R) D;
        auto q_lead_inv = q.lead().inv();
        while(R.deg() >= q.deg()) {
            D.a.push_back(R.lead() * q_lead_inv);
            if(D.lead() != 0) {
                for(size_t i = 1; i <= q.a.size(); i++) {
                    R.a[R.a.size() - i] -= D.lead() * q.a[q.a.size() - i];
                }
            }
            R.a.pop_back();
        }
        std::ranges::reverse(D.a);
        return std::pair{D, R};
    }
    auto divmod(auto const& p, auto const& q) {
        assert(!q.is_zero());
        int d = p.deg() - q.deg();
        auto D = decltype(p){};
        if(d >= 0) {
            D = (p.reverse().mod_xk(d + 1) * q.reverse().inv(d + 1)).mod_xk(d + 1).reverse(d + 1);
        }
        return std::pair{D, p - D * q};
    }
}
#endif // CP_ALGO_ALGEBRA_POLY_IMPL_DIV_HPP
