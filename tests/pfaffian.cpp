#include <bits/stdc++.h>
#include "cp-algo/linalg/matrix.hpp"
using namespace cp_algo::linalg;
using namespace cp_algo::math;
size_t checks = 0;
template<class T>
T naive(matrix<T> const& a, std::vector<size_t> const& ids) {
    if(ids.empty()) {
        return 1;
    }
    T ans = 0;
    for(size_t j = 1; j < ids.size(); j++) {
        std::vector<size_t> rest;
        for(size_t k = 1; k < ids.size(); k++) {
            if(k != j) {
                rest.push_back(ids[k]);
            }
        }
        T term = a[ids[0]][ids[j]] * naive(a, rest);
        ans += j % 2 ? term : -term;
    }
    return ans;
}
template<class T>
void check() {
    std::mt19937 rng(16288);
    for(size_t n = 0; n <= 10; n += 2) {
        for(int rep = 0; rep < 30; rep++) {
            matrix<T> a(n);
            for(size_t i = 0; i < n; i++) {
                for(size_t j = i + 1; j < n; j++) {
                    a[i][j] = rep % 3 == 0 ? rng() % 2 : rng() % 17;
                    a[j][i] = -a[i][j];
                }
            }
            if(rep % 5 == 0 && n) {
                for(size_t i = 0; i < n; i++) {
                    a[0][i] = a[i][0] = 0;
                }
            }
            std::vector<size_t> ids(n);
            std::iota(ids.begin(), ids.end(), 0);
            auto before = a;
            T p = a.pfaffian();
            assert(p == naive(a, ids));
            assert(p * p == a.det());
            assert(a == before);
            checks++;
        }
    }
}
template<class T>
void more() {
    // Regression: a pivot permutation must normalize lazy residues before negation.
    matrix<T> a(6);
    std::array<unsigned, 15> upper = {0, 0, 159749806, 963664571, 0,         0,         141540006, 0,
                                      0, 0, 689237731, 750435215, 669994770, 249911399, 332313983};
    size_t pos = 0;
    for(size_t i = 0; i < 6; i++) {
        for(size_t j = i + 1; j < 6; j++) {
            a[i][j] = upper[pos++];
            a[j][i] = -a[i][j];
        }
    }
    std::vector<size_t> ids(6);
    std::iota(ids.begin(), ids.end(), 0);
    assert(a.pfaffian() == naive(a, ids));
    checks++;
    std::mt19937 rng(7451);
    for(size_t n: {12, 16, 32, 64}) {
        for(int type = 0; type < 6; type++) {
            matrix<T> a(n);
            for(size_t i = 0; i < n; i++) {
                for(size_t j = i + 1; j < n; j++) {
                    a[i][j] = type == 0 ? T(-1) : (rng() % 4 == 0 ? T(rng()) : T(0));
                    a[j][i] = -a[i][j];
                }
            }
            if(type == 2) {
                for(size_t i = 0; i < n; i++) {
                    a[0][i] = a[i][0] = 0;
                }
            }
            if(type == 3) {
                for(size_t i = 0; i < n / 2; i++) {
                    for(size_t j = 0; j < n / 2; j++) {
                        a[i][j] = 0;
                    }
                }
            }
            auto p = a.pfaffian();
            assert(p * p == a.det());
            auto b = a;
            size_t u = 0, v = n - 1;
            std::swap(b[u], b[v]);
            for(size_t i = 0; i < n; i++) {
                std::swap(b[i][u], b[i][v]);
            }
            assert(b.pfaffian() == -p);
            checks++;
        }
    }
}
int main() {
    check<modint<998244353LL>>();
    check<modint<1000000007LL>>();
    check<modint<1073741789LL>>();
    check<modint<998244353>>();
    more<modint<998244353LL>>();
    more<modint<1073741789LL>>();
    more<modint<998244353>>();
    dynamic_modint<int64_t>::with_mod(998244353, [&] {
        check<dynamic_modint<int64_t>>();
        more<dynamic_modint<int64_t>>();
    });
    std::cout << checks << " Pfaffian cases passed signed matching expansion, singular and determinant checks\n";
}
