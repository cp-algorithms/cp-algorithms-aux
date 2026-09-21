#include "cp-algo/linalg/block_inverse.hpp"
#include <iostream>
#include <random>

using namespace cp_algo::linalg;
using namespace cp_algo::math;

template<typename base>
size_t check_inverse() {
    using M = matrix<base>;
    std::mt19937 rng(311);
    size_t checks = 0;
    for(size_t n: {127, 128, 129, 250, 257}) for(int type = 0; type < 6; type++) {
        M a(n);
        if(type <= 2) {
            for(auto &x: a.elements()) x = rng();
            if(type == 1) for(auto &x: a[0]) x = 0;
            if(type == 2) a[n - 1] = a[0];
        } else if(type == 3) {
            // Invertible matrix with a singular leading Schur block.
            for(size_t i = 0; i < n; i++) a[i][(i + n / 2) % n] = 1;
        } else if(type == 4) {
            for(size_t i = 0; i < n; i++) for(size_t j = 0; j < n; j++) {
                a[i][j] = (i + 1) * (j + 2);
            }
        } else {
            a = M::eye(n);
            a.gauss(); // Exercise input rows with cached pivots.
        }
        auto copy = a;
        auto [det, inverse] = block_inverse(a);
        auto [expected_det, expected_inverse] = a.inv();
        assert(a == copy && det == a.det());
        assert(det == expected_det && inverse == expected_inverse);
        if(det != base(0)) {
            assert(a * inverse == M::eye(n) && inverse * a == M::eye(n));
            assert(inverse.rank() == n); // Output rows must have valid pivot metadata.
        }
        checks++;
    }
    return checks;
}

int main() {
    auto checks = check_inverse<modint<998244353LL>>()
                + check_inverse<modint<998244353>>()
                + check_inverse<modint<1000000007LL>>();
    std::cout << checks << " Schur boundary, singular-leading-block, permutation and dense checks passed\n";
}
