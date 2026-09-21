#include "cp-algo/linalg/block_inverse.hpp"
#include <iostream>
#include <random>

using namespace cp_algo::linalg;
using namespace cp_algo::math;

size_t checks = 0;

template<typename base, typename row = modint_vec<base>>
void check_solutions() {
    using M = matrix<base, row>;
    std::mt19937 rng(821);
    for(size_t n: {0, 1, 2, 3, 7, 31, 32, 33, 65}) {
        for(size_t m: {0, 1, 2, 5, 31, 33, 65}) {
            // Empty matrices do not store a column count.
            if(n == 0 && m != 0) continue;
            for(int type = 0; type < 5; type++) {
                M a(n, m);
                for(auto &x: a.elements()) x = rng();
                if(type == 0) for(auto &x: a.elements()) x = 0;
                if(type == 1) for(size_t i = 1; i < n; i++) a[i] = a[0];
                if(type == 2) for(size_t i = 0; i < n; i++) {
                    for(size_t j = 0; j < m; j++) if(j % 3) a[i][j] = 0;
                }
                if(type == 3) for(size_t i = 0; i < n; i++) {
                    for(size_t j = 0; j < m; j++) a[i][j] = i + j + 1 == m;
                }
                auto kernel = a.kernel();
                size_t rank = a.rank();
                assert(kernel.n() == m - rank && kernel.rank() == kernel.n());
                if(kernel.n()) assert(a * kernel.T() == M(n, kernel.n()));
                for(size_t nrhs: {0, 1, 3}) {
                    if(n == 0 && nrhs != 0) continue;
                    M rhs(n, nrhs);
                    for(auto &x: rhs.elements()) x = type % 2 ? base(0) : base(rng());
                    auto result = a.solve(rhs);
                    assert(bool(result) == ((a | rhs).rank() == rank));
                    if(result) {
                        auto const& [solution, basis] = *result;
                        assert(solution.n() == nrhs);
                        // Scalar residuals also cover systems with zero variables.
                        for(size_t i = 0; i < n; i++) for(size_t j = 0; j < nrhs; j++) {
                            base value = 0;
                            for(size_t k = 0; k < m; k++) value += a[i][k] * solution[j][k];
                            assert(value == rhs[i][j]);
                        }
                        assert(basis.n() == m - rank && basis.rank() == basis.n());
                        if(basis.n()) assert(a * basis.T() == M(n, basis.n()));
                    }
                    checks++;
                }
                if(n == m) {
                    auto [det, inverse] = a.inv();
                    assert(det == a.det());
                    assert((det != base(0)) == (rank == n));
                    if(det != base(0)) {
                        assert(a * inverse == M::eye(n) && inverse * a == M::eye(n));
                    }
                    auto [block_det, block_inv] = block_inverse(a);
                    assert(block_det == det && block_inv == inverse);
                    checks++;
                }
            }
        }
    }
}

int main() {
    check_solutions<modint<998244353LL>>();
    check_solutions<modint<998244353>>();
    check_solutions<modint<1000000007LL>>();
    check_solutions<modint<17>>();
    check_solutions<modint<2>>();
    check_solutions<modint<998244353LL>, vec<modint<998244353LL>>>();
    dynamic_modint<int64_t>::with_mod(998244353, [] {
        check_solutions<dynamic_modint<int64_t>>();
    });
    std::cout << checks << " kernel, multiple-RHS solve and inverse property checks passed\n";
}
