#include "cp-algo/math/newton.hpp"
#include "cp-algo/math/poly/compose.hpp"
#include <cassert>
#include <cstddef>
#include <limits>
#include <random>
#include <utility>
using namespace cp_algo::math;
using T = modint<998244353>;
using P = poly_t<T>;

int main() {
    std::mt19937 rng(42);
    for(size_t n: {0, 1, 2, 3, 31, 32, 33, 63, 64, 65, 129}) {
        P::Vector a(n);
        for(auto& x: a) {x = 1 + rng() % 100;}
        P p(a);
        auto q = newton([&](P const& f, size_t m) {
            return std::pair{(f * p - P(1)).mod_xk(m), p.mod_xk(m)};
        }, n ? a[0].inv() : T(1), n);
        P::Vector expected(n);
        if(n) {expected[0] = a[0].inv();}
        for(size_t i = 1; i < n; ++i) {
            for(size_t j = 1; j <= i; ++j) {expected[i] -= a[j] * expected[i - j];}
            expected[i] *= expected[0];
        }
        assert(q == P(expected));
        auto original = p;
        assert(p.bisect(std::numeric_limits<size_t>::max()) == p.bisect());
        p += p;
        assert(p == original * T(2));
        p -= p;
        assert(p.is_zero());
        assert(compose_large(original, P(P::Vector{0, 0, 1}), n) == original.x2().mod_xk(n));
        assert(compose_large(P(), original, n).is_zero());
    }
}
