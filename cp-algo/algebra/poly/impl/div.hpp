#ifndef CP_ALGO_ALGEBRA_POLY_IMPL_DIV_HPP
#define CP_ALGO_ALGEBRA_POLY_IMPL_DIV_HPP
#include "../../common.hpp"
#include <cassert>
// operations related to polynomial division
namespace cp_algo::algebra::poly::impl {
    auto divmod_slow(auto const& p, auto const& q) {
        auto R = p;
        auto D = decltype(p){};
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
        R.normalize();
        return std::array{D, R};
    }
    template<typename poly>
    auto divmod_hint(poly const& p, poly const& q, poly const& qri) {
        assert(!q.is_zero());
        int d = p.deg() - q.deg();
        if(std::min(d, q.deg()) < magic) {
            return divmod_slow(p, q);
        }
        poly D;
        if(d >= 0) {
            D = (p.reverse().mod_xk(d + 1) * qri.mod_xk(d + 1)).mod_xk(d + 1).reverse(d + 1);
        }
        return std::array{D, p - D * q};
    }
    auto divmod(auto const& p, auto const& q) {
        assert(!q.is_zero());
        int d = p.deg() - q.deg();
        if(std::min(d, q.deg()) < magic) {
            return divmod_slow(p, q);
        }
        return divmod_hint(p, q, q.reverse().inv(d + 1));
    }

    template<typename poly>
    poly powmod_hint(poly const& p, int64_t k, poly const& md, poly const& mdri) {
        return bpow(p, k, poly(1), [&](auto const& p, auto const& q){
            return divmod_hint(p * q, md, mdri)[1];
        });
    }
    template<typename poly>
    auto powmod(poly const& p, uint64_t k, poly const& md) {
        int d = md.deg();
        if(md == poly::xk(d)) {
            return p.pow(k, d);
        }
        if(md == poly::xk(d) - poly(1)) {
            return p.powmod_circular(k, d);
        }
        auto mdinv = md.reverse().inv(md.deg() + 1);
        return powmod_hint(p, k, md, mdinv);
    }
}
#endif // CP_ALGO_ALGEBRA_POLY_IMPL_DIV_HPP
