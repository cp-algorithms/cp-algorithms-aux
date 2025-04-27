#ifndef CP_ALGO_MATH_POLY_IMPL_DIV_HPP
#define CP_ALGO_MATH_POLY_IMPL_DIV_HPP
#include "../../fft.hpp"
#include "../../common.hpp"
#include <cassert>
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
            D = (p.reversed().mod_xk(d + 1) * qri.mod_xk(d + 1)).mod_xk(d + 1).reversed(d + 1);
        }
        return std::array{D, p - D * q};
    }
    auto divmod(auto const& p, auto const& q) {
        assert(!q.is_zero());
        int d = p.deg() - q.deg();
        if(std::min(d, q.deg()) < magic) {
            return divmod_slow(p, q);
        }
        return divmod_hint(p, q, q.reversed().inv(d + 1));
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
                auto mdr = md.reversed();
                return (mdr.inv(k - md.deg() + 1, md.deg()) * mdr).reversed(md.deg());
            }
        }
        if(md == poly::xk(d)) {
            return p.pow(k, d);
        }
        if(md == poly::xk(d) - poly(1)) {
            return p.powmod_circular(k, d);
        }
        return powmod_hint(p, k, md, md.reversed().inv(md.deg() + 1));
    }
    template<typename poly>
    poly& inv_inplace(poly& q, int64_t k, size_t n) {
        using poly_t = std::decay_t<poly>;
        using base = poly_t::base;
        if(k <= std::max<int64_t>(n, size(q.a))) {
            return q.inv_inplace(k + n).div_xk_inplace(k);
        }
        if(k % 2) {
            return inv_inplace(q, k - 1, n + 1).div_xk_inplace(1);
        }
        auto [q0, q1] = q.bisect();
        auto qq = q0 * q0 - (q1 * q1).mul_xk_inplace(1);
        inv_inplace(qq, k / 2 - q.deg() / 2, (n + 1) / 2 + q.deg() / 2);
        size_t N = fft::com_size(size(q0.a), size(qq.a));
        auto q0f = fft::dft<base>(q0.a, N);
        auto q1f = fft::dft<base>(q1.a, N);
        auto qqf = fft::dft<base>(qq.a, N);
        size_t M = q0.deg() + (n + 1) / 2;
        typename poly::Vector A, B;
        q0f.mul(qqf, A, M);
        q1f.mul_inplace(qqf, B, M);
        q.a.resize(n + 1);
        for(size_t i = 0; i < n; i += 2) {
            q.a[i] = A[q0.deg() + i / 2];
            q.a[i + 1] = -B[q0.deg() + i / 2];
        }
        q.a.pop_back();
        q.normalize();
        return q;
    }
    template<typename poly>
    poly& inv_inplace(poly& p, size_t n) {
        using poly_t = std::decay_t<poly>;
        using base = poly_t::base;
        if(n == 1) {
            return p = base(1) / p[0];
        }
        // Q(-x) = P0(x^2) + xP1(x^2)
        auto [q0, q1] = p.bisect(n);
        
        size_t N = fft::com_size(size(q0.a), (n + 1) / 2);
        
        auto q0f = fft::dft<base>(q0.a, N);
        auto q1f = fft::dft<base>(q1.a, N);

        // Q(x)*Q(-x) = Q0(x^2)^2 - x^2 Q1(x^2)^2
        auto qq = poly_t(q0f * q0f) - poly_t(q1f * q1f).mul_xk_inplace(1);

        inv_inplace(qq, (n + 1) / 2);
        auto qqf = fft::dft<base>(qq.a, N);
        
        typename poly::Vector A, B;
        q0f.mul(qqf, A, (n + 1) / 2);
        q1f.mul_inplace(qqf, B, (n + 1) / 2);
        p.a.resize(n + 1);
        for(size_t i = 0; i < n; i += 2) {
            p.a[i] = A[i / 2];
            p.a[i + 1] = -B[i / 2];
        }
        p.a.pop_back();
        p.normalize();
        return p;
    }
}
#endif // CP_ALGO_MATH_POLY_IMPL_DIV_HPP
