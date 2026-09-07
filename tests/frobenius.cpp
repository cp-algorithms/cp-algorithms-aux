#include "cp-algo/linalg/frobenius.hpp"
#include <iostream>

using namespace cp_algo::math;
using namespace cp_algo::linalg;

template<typename base>
size_t check_decomposition() {
    using M = matrix<base>;
    size_t cases = 0;
    for(size_t n: {0, 1, 2, 3, 7, 16, 17, 33}) {
        for(int type = 0; type < 6; type++) {
            M a(n);
            if(type == 0) a = M::random(n);
            if(type == 1) a = M::eye(n);
            if(type == 2) for(size_t i = 0; i < n; i++) a[i][i] = i % 3;
            if(type == 3) for(size_t i = 0; i + 1 < n; i++) a[i][i + 1] = 1;
            if(type == 4 && n) a = M::random(n, 2) * M::random(2, n);
            auto [T, Ti, blocks] = frobenius_form<full>(a);
            M C(n);
            size_t start = 0;
            for(auto const& p: blocks) {
                size_t d = p.deg();
                assert(d > 0 && start + d <= n);
                for(size_t j = 0; j + 1 < d; j++) C[start + j][start + j + 1] = 1;
                for(size_t j = 0; j < d; j++) C[start + d - 1][start + j] = -p[j] / p[d];
                start += d;
            }
            assert(start == n && Ti * T == M::eye(n));
            assert(T * a == C * T);
            auto only_blocks = frobenius_form(a);
            for(int x: {0, 1, 2, -1}) {
                base expected = (M::eye(n) * base(x) - a).det();
                base value = 1, full_value = 1;
                for(auto const& p: only_blocks) value *= p.eval(base(x));
                for(auto const& p: blocks) full_value *= p.eval(base(x));
                assert(value == expected && full_value == expected);
            }
            cases++;
        }
    }
    return cases;
}

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
    auto decompositions = check_decomposition<base>() + check_decomposition<modint<1000000007LL>>();
    std::cout << decompositions << " Frobenius decompositions and characteristic-polynomial checks passed\n";
}
