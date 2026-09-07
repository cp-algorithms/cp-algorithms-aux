#include <bits/stdc++.h>
#include "cp-algo/linalg/matrix.hpp"
using namespace cp_algo::linalg;
using namespace cp_algo::math;
template<class T> T naive(matrix<T> const& a, std::vector<size_t> const& ids) {
    if(ids.empty()) return 1;
    T ans = 0;
    for(size_t j = 1; j < ids.size(); j++) {
        std::vector<size_t> rest;
        for(size_t k = 1; k < ids.size(); k++) if(k != j) rest.push_back(ids[k]);
        T term = a[ids[0]][ids[j]] * naive(a, rest);
        ans += j % 2 ? term : -term;
    }
    return ans;
}
template<class T> void check() {
    std::mt19937 rng(16288);
    for(size_t n = 0; n <= 10; n += 2) for(int rep = 0; rep < 30; rep++) {
        matrix<T> a(n);
        for(size_t i = 0; i < n; i++) for(size_t j = i + 1; j < n; j++) {
            a[i][j] = rep % 3 == 0 ? rng() % 2 : rng() % 17;
            a[j][i] = -a[i][j];
        }
        if(rep % 5 == 0 && n) for(size_t i = 0; i < n; i++) a[0][i] = a[i][0] = 0;
        std::vector<size_t> ids(n);
        std::iota(ids.begin(), ids.end(), 0);
        auto before = a;
        T p = a.pfaffian();
        assert(p == naive(a, ids));
        assert(p * p == a.det());
        assert(a == before);
    }
}
int main() {
    check<modint<998244353LL>>();
    check<modint<1000000007LL>>();
    std::cout << "Pfaffians passed signed matching expansion, singular and determinant checks\n";
}
