#include "cp-algo/math/poly/recurrence.hpp"
#include <cassert>
#include <random>
using namespace cp_algo::math;
using T=modint<998244353>; using P=poly_t<T>;
int main() {
    assert(kth_rec(P::xk(100), P(1), 100)==T(1));
    assert(kth_rec(P::xk(100), P(1), 1000000000000000000LL)==T(0));
    std::mt19937 rng(451);
    for(int a=1;a<=40;a++) for(int b=1;b<=10;b++) {
        P::Vector av(a),bv(b);
        for(auto &x:av)x=rng()%100;
        for(auto &x:bv)x=rng()%100;
        bv[0]=1;
        P A(av),B(bv); auto want=(A*inv(B,100)).mod_xk(100);
        for(int k:{0,1,7,30,70,99})assert(kth_rec(A,B,k)==want[k]);
    }
}
