#include "cp-algo/linalg/frobenius.hpp"
#include <iostream>

using namespace cp_algo::math;
using namespace cp_algo::linalg;

int main() {
    using base = modint<998244353LL>;
    using M = matrix<base>;
    cp_algo::random::gen.seed(613);
    size_t cases = 0;
    for(size_t n: {0, 1, 3, 17, 64}) {
        for(int type = 0; type < 3; type++) {
            M a = type == 0 ? M::random(n) : type == 1 ? M::eye(n) : M(n);
            M expected = M::eye(n);
            for(uint64_t k = 0; k <= 3; k++) {
                assert(frobenius_pow(a, k) == expected);
                cases++;
                M next(n);
                for(size_t i = 0; i < n; i++)
                for(size_t j = 0; j < n; j++)
                for(size_t t = 0; t < n; t++) next[i][j] += expected[i][t] * a[t][j];
                expected = std::move(next);
            }
        }
    }
    std::cout << cases << " small Frobenius powers passed against scalar multiplication\n";
}
