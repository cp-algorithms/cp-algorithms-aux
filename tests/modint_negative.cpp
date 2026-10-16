#include "cp-algo/number_theory/modint.hpp"
#include <cassert>
#include <cstdint>
#include <limits>
using namespace cp_algo::math;
template <class T> void check() {
    for(int64_t x : {int64_t(-100), int64_t(-1000000000000000000LL), std::numeric_limits<int64_t>::min(), int64_t(-1),
                     int64_t(0)}) {
        __int128 r = __int128(x) % T::mod();
        if(r < 0)
            r += T::mod();
        assert(T(x).getr() == r);
    }
}
int main() {
    check<modint<7>>();
    check<modint<998244353>>();
    check<modint<2305843009213693951LL>>();
    dynamic_modint<int>::with_mod(7, [] { check<dynamic_modint<int>>(); });
    dynamic_modint<int64_t>::with_mod(998244353, [] { check<dynamic_modint<int64_t>>(); });
}
