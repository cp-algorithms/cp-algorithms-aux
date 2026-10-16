#include "cp-algo/math/bigint.hpp"
#include <cassert>
#include <limits>
#include <sstream>
using namespace cp_algo::math;
template <base_v B = x10> bigint<B> parse(std::string text) {
    return bigint<B>(std::span<char>(text));
}
template <base_v B> std::string str(bigint<B> const &x) {
    std::ostringstream s;
    print_bigint(static_cast<std::ostream &>(s), x);
    return s.str();
}
int main() {
    for(int64_t x : {int64_t(0), int64_t(1), int64_t(-1), int64_t(10000000000000000LL), int64_t(-10000000000000000LL),
                     std::numeric_limits<int64_t>::min(), std::numeric_limits<int64_t>::max()}) {
        auto want = std::to_string(x);
        assert(str(bigint<>(x)) == want);
        assert((bigint<>(x) <=> parse(want)) == 0);
    }
    assert(str(bigint<x16>(std::numeric_limits<int64_t>::min())) == "-8000000000000000");
    assert(str(bigint<x16>(std::numeric_limits<int64_t>::max())) == "7FFFFFFFFFFFFFFF");
    assert(str(parse("9999999999999999") * int64_t(1000000000000000000LL)) == "9999999999999999000000000000000000");
    assert(str(parse("-9999999999999999") * int64_t(1000000000000000000LL)) == "-9999999999999999000000000000000000");
    assert(str(bigint<>(2) * std::numeric_limits<int64_t>::min()) == "-18446744073709551616");
    assert(str(bigint<>(0) * int64_t(-1)) == "0");

    std::optional<bigint<>> default_value(std::in_place, -1);
    default_value.reset();
    default_value.emplace();
    assert(str(*default_value) == "0");
    assert((*default_value <=> bigint<>(0)) == 0);
    *default_value += bigint<>(2);
    assert(str(*default_value) == "2");
}
