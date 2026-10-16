#include "cp-algo/number_theory/modint.hpp"
#include <cassert>
#include <sstream>
#include <vector>

using namespace cp_algo::math;

template <class T> void check() {
    std::vector<typename T::UInt> values{0, 1, typename T::UInt(T::mod() / 2), typename T::UInt(T::mod() - 1)};
    std::stringstream input;
    for(auto x : values)
        input << x << ' ';
    for(auto x : values) {
        T value;
        input >> value;
        assert(value.getr() == x);
        assert((value + T(1)).getr() == (x + 1) % T::mod());
    }
}

int main() {
    check<modint<998244353>>();
    check<modint<2305843009213693951LL>>();
    for(int m : {2, 998244353, 2147483647}) {
        dynamic_modint<>::with_mod(m, [] { check<dynamic_modint<>>(); });
    }
    for(int64_t m : {int64_t(998244353), int64_t(2305843009213693951)}) {
        dynamic_modint<int64_t>::with_mod(m, [] { check<dynamic_modint<int64_t>>(); });
    }
}
