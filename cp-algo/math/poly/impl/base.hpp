#ifndef CP_ALGO_MATH_POLY_IMPL_BASE_HPP
#define CP_ALGO_MATH_POLY_IMPL_BASE_HPP
#include <functional>
#include <algorithm>
#include <iostream>
// really basic operations, typically taking O(n)
namespace cp_algo::math::poly::impl {
    template<typename polyn>
    void normalize(polyn& p) {
        while(p.deg() >= 0 && p.lead() == typename polyn::base(0)) {
            p.a.pop_back();
        }
    }
    auto neg_inplace(auto &&p) {
        std::ranges::transform(p.a, begin(p.a), std::negate{});
        return p;
    }
    auto& scale(auto &p, auto x) {
        for(auto &it: p.a) {
            it *= x;
        }
        p.normalize();
        return p;
    }
    auto& add(auto &p, auto const& q) {
        p.a.resize(std::max(p.a.size(), q.a.size()));
        std::ranges::transform(p.a, q.a, begin(p.a), std::plus{});
        normalize(p);
        return p;
    }
    auto& sub(auto &p, auto const& q) {
        p.a.resize(std::max(p.a.size(), q.a.size()));
        std::ranges::transform(p.a, q.a, begin(p.a), std::minus{});
        normalize(p);
        return p;
    }
    auto substr(auto const& p, size_t l, size_t k) {
        return std::vector(
            begin(p.a) + std::min(l, p.a.size()),
            begin(p.a) + std::min(l + k, p.a.size())
        );
    }
    auto& reverse(auto &p, size_t n) {
        p.a.resize(n);
        std::ranges::reverse(p.a);
        normalize(p);
        return p;
    }
}
#endif // CP_ALGO_MATH_POLY_IMPL_BASE_HPP
