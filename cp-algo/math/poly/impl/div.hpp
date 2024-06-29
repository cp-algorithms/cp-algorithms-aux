#ifndef CP_ALGO_MATH_POLY_IMPL_DIV_HPP
#define CP_ALGO_MATH_POLY_IMPL_DIV_HPP
#include "../../fft.hpp"
#include "../../common.hpp"
#include <cassert>
#include <deque>
// operations related to polynomial division
namespace cp_algo::math::poly::impl {
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
    auto powmod(poly const& p, int64_t k, poly const& md) {
        int d = md.deg();
        if(p == poly::xk(1) && false) { // does it actually speed anything up?..
            if(k < md.deg()) {
                return poly::xk(k);
            } else {
                auto mdr = md.reverse();
                return (mdr.inv(k - md.deg() + 1, md.deg()) * mdr).reverse(md.deg());
            }
        }
        if(md == poly::xk(d)) {
            return p.pow(k, d);
        }
        if(md == poly::xk(d) - poly(1)) {
            return p.powmod_circular(k, d);
        }
        return powmod_hint(p, k, md, md.reverse().inv(md.deg() + 1));
    }

    auto interleave(auto const& p) {
        auto [p0, p1] = p.bisect();
        return p0 * p0 - (p1 * p1).mul_xk(1);
    }
    template<typename poly>
    poly inv(poly const& q, int64_t k, size_t n) {
        if(k <= std::max<int64_t>(n, size(q.a))) {
            return q.inv(k + n).div_xk(k);
        }
        if(k % 2) {
            return inv(q, k - 1, n + 1).div_xk(1);
        }
        
        auto qq = inv(interleave(q), k / 2 - q.deg() / 2, (n + 1) / 2 + q.deg() / 2);
        auto [q0, q1] = q.negx().bisect();
        return (
            (q0 * qq).x2() + (q1 * qq).x2().mul_xk(1)
        ).div_xk(2*q0.deg()).mod_xk(n);
    }
    template<typename poly>
    void inv_inplace(poly&& p, size_t n) {
        using poly_t = std::decay_t<poly>;
        using base = poly_t::base;
        if(n == 1) {
            p = base(1) / p[0];
            return;
        }
        // Q(-x) = P0(x^2) + xP1(x^2)
        auto [q0, q1] = p.bisect(n);
        
        int N = fft::com_size((n + 1) / 2, (n + 1) / 2);
        
        auto q0f = fft::dft<base>(q0.a, N);
        auto q1f = fft::dft<base>(q1.a, N);

        // Q(x)*Q(-x) = Q0(x^2)^2 - x^2 Q1(x^2)^2
        auto qq = poly_t(q0f * q0f) - poly_t(q1f * q1f).mul_xk(1);
        inv_inplace(qq, (n + 1) / 2);
        auto qqf = fft::dft<base>(qq.a, N);
        
        auto A = q0f * qqf;
        auto B = q1f * qqf;
        p.a.resize(n + 1);
        for(size_t i = 0; i < n; i += 2) {
            p.a[i] = A[i / 2];
            p.a[i + 1] = -B[i / 2];
        }
        p.a.pop_back();
    }
}
#endif // CP_ALGO_MATH_POLY_IMPL_DIV_HPP
