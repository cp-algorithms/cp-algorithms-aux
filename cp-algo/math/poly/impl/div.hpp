#ifndef CP_ALGO_MATH_POLY_IMPL_DIV_HPP
#define CP_ALGO_MATH_POLY_IMPL_DIV_HPP
#include "../inv.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::poly::impl {
    template<typename T>
    std::array<poly_t<T>, 2> divmod_slow(poly_t<T> p, poly_t<T> const& q) {
        poly_t<T> d;
        auto qi = q.lead().inv();
        while(p.deg() >= q.deg()) {
            d.a.push_back(p.lead() * qi);
            if(d.lead() != T(0)) {
                for(size_t i = 1; i <= q.a.size(); i++) {
                    p.a[p.a.size() - i] -= d.lead() * q.a[q.a.size() - i];
                }
            }
            p.a.pop_back();
        }
        std::ranges::reverse(d.a);
        p.normalize();
        return {std::move(d), std::move(p)};
    }
    template<typename T>
    std::array<poly_t<T>, 2> divmod_hint(poly_t<T> p, poly_t<T> const& q, poly_t<T> const& qri) {
        assert(!q.is_zero());
        int n = p.deg() - q.deg();
        if(std::min(n, q.deg()) < magic) {
            return divmod_slow(std::move(p), q);
        }
        auto d = p.reversed().mod_xk(n + 1);
        d.mul_truncate(qri.mod_xk(n + 1), n + 1).reverse(n + 1);
        p -= d * q;
        return {std::move(d), std::move(p)};
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_POLY_IMPL_DIV_HPP
