#ifndef CP_ALGO_ALGEBRA_POLY_IMPL_BASE_HPP
#define CP_ALGO_ALGEBRA_POLY_IMPL_BASE_HPP
#include <functional>
#include <algorithm>
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
}
#endif // CP_ALGO_ALGEBRA_POLY_IMPL_BASE_HPP
