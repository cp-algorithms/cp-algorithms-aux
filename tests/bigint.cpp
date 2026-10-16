#include "cp-algo/math/bigint.hpp"
#include <random>
using namespace cp_algo;
using namespace cp_algo::math;

template<base_v radix>
void check() {
    using B = bigint<radix>;
    std::mt19937_64 rng(71823);
    auto value = [](B const& a) {
        __int128 result = 0;
        for(auto digit: a.digits | std::views::reverse) {result = result * B::Base + digit;}
        return a.negative ? -result : result;
    };
    std::vector<B> cases{B(0), B(1), B(-1), B({0, 1}, false), B({0, 1}, true)};
    auto zero = B(0);
    zero.negate(); // A negative zero still compares equal to zero.
    cases.push_back(zero);
    for(int i = 0; i < 100; i++) {
        cases.emplace_back(big_basic_string<uint64_t>{rng() % B::Base, rng() % B::Base}, rng() % 2);
    }
    for(auto const& a: cases) {
        for(auto const& b: cases) {
            assert((a <=> b) == (value(a) <=> value(b)));
        }
        for(size_t n = 0; n <= 3; n++) {
            __int128 power = 1;
            // Two limbs fit in __int128; larger truncations leave the input unchanged.
            for(size_t i = 0; i < std::min(n, size_t(2)); i++) {power *= B::Base;}
            assert(value(a.take(n)) == (n >= 2 ? value(a) : value(a) % power));
            assert(value(a.drop(n)) == (n >= 2 ? 0 : value(a) / power));
        }
    }
    for(int64_t x: {-1000000, -1, 0, 1, 1000000}) {
        B a(x);
        assert(value(a.pad(1)) == __int128(x) * B::Base);
        assert(a.pad(1).drop(1).digits == a.digits);
    }
}
int main() {
    check<x10>();
    check<x16>();
    std::cout << "Bigint ordering and slicing passed independent integer checks\n";
}
