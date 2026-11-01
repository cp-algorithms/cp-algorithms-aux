#include "cp-algo/math/poly/transform.hpp"
#include <cassert>
#include <random>
using namespace cp_algo::math;
template<class T> void check() {
    std::mt19937 rng(42);
    for(int n = 0; n <= 48; ++n) {
        typename poly_t<T>::Vector a(n);
        for(auto &x : a) x = rng() % T::mod();
        poly_t<T> p(a);
        auto q = prefix_sum(p);
        assert(q[0] == T(0));
        for(int x = -5; x <= 55; ++x) assert(q.eval(T(x + 1)) - q.eval(T(x)) == p.eval(T(x)));
    }
    assert(prefix_sum(poly_t<T>{}).is_zero());
}
int main() {check<modint<998244353>>(); check<modint<101>>();}
