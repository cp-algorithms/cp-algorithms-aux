#include "cp-algo/math/decimal.hpp"
#include <cassert>
#include <limits>
#include <sstream>
using namespace cp_algo::math;
bigint<> parse(std::string text) {
    return bigint<>(std::span<char>(text));
}
std::string str(bigint<> const &x) {
    std::ostringstream s;
    print_bigint(static_cast<std::ostream &>(s), x);
    return s.str();
}
int main() {
    for(int a = -30; a <= 30; a++)
        for(int b = -9; b <= 9; b++)
            if(b) {
                auto [q, r] = divmod(bigint<>(a), bigint<>(b));
                assert(str(q) == std::to_string(a / b) && str(r) == std::to_string(a % b));
            }
    for(bool negative_a : {false, true})
        for(bool negative_b : {false, true}) {
            auto a = parse("1234567890123456789012345678901234567890");
            auto b = parse("987654321098765432109");
            if(negative_a)
                a.negate();
            if(negative_b)
                b.negate();
            auto [q, r] = divmod(a, b);
            assert((q * b + r <=> a) == 0);
            auto ar = r, ab = b;
            ar.negative = ab.negative = false;
            assert(ar < ab && (r.digits.empty() || r.negative == negative_a));
        }
    auto [q, r] = divmod_fast(parse("18446744073709551617"), std::numeric_limits<int64_t>::min());
    assert(str(q) == "-2" && str(r) == "1");
    for(int q = -40; q <= 40; q++) {
        decimal<> x(int64_t(q) * 2500000000000000LL, -1);
        auto want = q < 0 ? -((-q + 2) / 4) : (q + 2) / 4;
        assert(str(x.round()) == std::to_string(want));
    }

    for(int scale = -3; scale <= 3; scale++)
        for(int v : {-7, -2, 1, 2, 7}) {
            decimal<> x(v, scale);
            auto reciprocal = x.inv(3);
            auto error = decimal<>(1) - x * reciprocal;
            assert(error.magnitude() <= -3);
        }

    auto printed = [](decimal<> x) {
        std::ostringstream out;
        auto previous = std::cout.rdbuf(out.rdbuf());
        x.print();
        std::cout.rdbuf(previous);
        return out.str();
    };
    assert(printed(decimal<>(1)) == "1\n");
    assert(printed(decimal<>(0, -1)) == "0\n");
    assert(printed(decimal<>(-7500000000000000LL, -1)) == "-0.7500000000000000\n");
    assert(printed(decimal<>(-12500000000000000LL, -1)) == "-1.2500000000000000\n");
    assert(printed(decimal<>(1, -2)) == "0." + std::string(31, '0') + "1\n");
    assert(printed(decimal<>(1, 1)) == "10000000000000000\n");
}
