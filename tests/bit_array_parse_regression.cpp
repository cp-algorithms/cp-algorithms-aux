#include <cassert>
#include <bits/stdc++.h>
#ifdef TEST_MINIFIED
#include "cp-algo/min/structures/bit_array_util.hpp"
#else
#include "cp-algo/structures/bit_array_util.hpp"
#endif
template <class Bits> void check() {
    Bits bits;
    auto parse = [&](std::string input) {
        cp_algo::structures::from_string(bits, std::span<char>(input));
        auto actual = cp_algo::structures::to_string(bits);
        assert(std::string(actual.begin(), actual.end()) == input);
        assert(cp_algo::structures::count(bits) == size_t(std::count(input.begin(), input.end(), '1')));
    };
    parse("1");
    parse("0");
    parse("");
    for(size_t n : {1, 2, 63, 64, 65, 127, 128, 129, 255, 256, 257}) {
        parse(std::string(n, '1'));
        parse(std::string(n, '0'));
    }
    std::mt19937 rng(4261);
    for(int test = 0; test < 1000; ++test) {
        std::string input(1 + rng() % 257, '0');
        for(char &c : input)
            c += rng() % 2;
        parse(input);
    }
}
int main() {
    check<cp_algo::structures::dynamic_bit_array>();
    check<cp_algo::structures::bit_array<320>>();
}
