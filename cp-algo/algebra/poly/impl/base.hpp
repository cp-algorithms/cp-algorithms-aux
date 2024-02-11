#ifndef CP_ALGO_ALGEBRA_POLY_IMPL_BASE_HPP
#define CP_ALGO_ALGEBRA_POLY_IMPL_BASE_HPP
#include <functional>
#include <algorithm>
#include <iostream>
// really basic operations, typically taking O(n)
namespace cp_algo::algebra::poly::impl {
    void normalize(auto& p) {
        while(p.deg() >= 0 && p.lead() == 0) {
            p.a.pop_back();
        }
    }
    auto neg(auto p) {
        std::ranges::transform(p.a, begin(p.a), std::negate<>{});
        return p;
    }
    auto& add(auto &p, auto q) {
        p.a.resize(std::max(p.a.size(), q.a.size()));
        std::ranges::transform(p.a, q.a, begin(p.a), std::plus<>{});
        normalize(p);
        return p;
    }
    auto& sub(auto &p, auto q) {
        p.a.resize(std::max(p.a.size(), q.a.size()));
        std::ranges::transform(p.a, q.a, begin(p.a), std::minus<>{});
        normalize(p);
        return p;
    }
    auto mod_xk(auto const& p, size_t k) {
        return std::vector(begin(p.a), begin(p.a) + std::min(k, p.a.size()));
    }
    auto mul_xk(auto p, size_t k) {
        p.a.insert(begin(p.a), k, 0);
        normalize(p);
        return p;
    }
    auto div_xk(auto const& p, size_t k) {
        return std::vector(begin(p.a) + std::min(k, p.a.size()), end(p.a));
    }
    auto substr(auto const& p, size_t l, size_t r) {
        return std::vector(
            begin(p.a) + std::min(l, p.a.size()),
            begin(p.a) + std::min(r, p.a.size())
        );
    }
    auto reverse(auto p, size_t n) {
        p.a.resize(n);
        std::ranges::reverse(p.a);
        normalize(p);
        return p;
    }
}
#endif // CP_ALGO_ALGEBRA_POLY_IMPL_BASE_HPP
