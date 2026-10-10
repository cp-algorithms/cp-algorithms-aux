#include "cp-algo/linalg/matrix.hpp"
#include <random>
#include <iostream>
using namespace cp_algo;
using namespace cp_algo::math;
using namespace cp_algo::linalg;
template<auto mod> size_t check() {
    using B = modint<mod>;
    using M = matrix<B>;
    using P = cp_algo::linalg::impl::strassen_product<mod>;
    std::mt19937 rng(820471);
    std::array<size_t, 3> shapes[] = {{1,1,1},{3,7,9},{31,32,33},{33,65,129},
        {63,64,65},{64,64,64},{65,65,65},{97,129,33},{128,128,128},
        {129,129,129},{257,65,193},{257,1,64}};
    size_t cases = 0;
    for(auto [n,m,k]: shapes) for(int type = 0; type < 3; type++) for(int sign: {-1,0,1}) {
        M a(n,m), b(m,k), prior(n,k);
        for(auto &v: a.elements()) v = type == 0 ? B(rng()) : type == 1 ? B(0) : B(-1);
        for(auto &v: b.elements()) v = type == 0 ? B(rng()) : type == 1 ? B(0) : B(-1);
        for(auto &v: prior.elements()) v = type == 0 ? B(rng()) : type == 1 ? B(-1) : B(0);
        auto target = prior;
        auto c = sign == 1 ? P::template product<1>(a,b,std::move(target)) :
                 sign == -1 ? P::template product<-1>(a,b,std::move(target)) : P::product(a,b);
        assert(c.n() == n && c.m() == k);
        B scale = sign == -1 ? -1 : 1;
        if(type || std::max({n,m,k}) <= 65 || m == 1) {
            for(size_t i = 0; i < n; i++) for(size_t j = 0; j < k; j++) {
                B value = 0;
                if(type == 2) value = m;
                else if(type == 0) for(size_t t = 0; t < m; t++) value += a[i][t] * b[t][j];
                value *= scale;
                if(sign) value += prior[i][j];
                assert(c[i][j] == value);
            }
        } else for(int rep = 0; rep < 3; rep++) {
            std::vector<B> x(k), bx(m), expected(n), actual(n);
            for(auto &v: x) v = rng();
            for(size_t i = 0; i < m; i++) for(size_t j = 0; j < k; j++) bx[i] += b[i][j] * x[j];
            for(size_t i = 0; i < n; i++) {
                for(size_t j = 0; j < m; j++) expected[i] += a[i][j] * bx[j];
                expected[i] *= scale;
                for(size_t j = 0; j < k; j++) {
                    if(sign) expected[i] += prior[i][j] * x[j];
                    actual[i] += c[i][j] * x[j];
                }
            }
            assert(expected == actual);
        }
        cases++;
    }
    return cases;
}
int main() {
    auto count = check<998244353>() + check<1000000007>() + check<1073741789>();
    std::cout << count << " product/accumulate/subtract reference checks passed\n";
}
