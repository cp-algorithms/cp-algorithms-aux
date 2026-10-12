#include "cp-algo/math/bigint.hpp"
#include <cassert>
#include <limits>
#include <sstream>
using namespace cp_algo::math;
template<base_v B = x10> bigint<B> parse(std::string text) {return bigint<B>(std::span<char>(text));}
template<base_v B> std::string str(bigint<B> const& x) {std::ostringstream s; print_bigint(static_cast<std::ostream&>(s), x); return s.str();}
int main() {
    for(int64_t x: {int64_t(0), int64_t(1), int64_t(-1), int64_t(10000000000000000LL), int64_t(-10000000000000000LL), std::numeric_limits<int64_t>::min(), std::numeric_limits<int64_t>::max()}) {
        auto want = std::to_string(x);
        assert(str(bigint<>(x)) == want);
        assert((bigint<>(x) <=> parse(want)) == 0);
    }
    assert(str(bigint<x16>(std::numeric_limits<int64_t>::min())) == "-8000000000000000");
    assert(str(bigint<x16>(std::numeric_limits<int64_t>::max())) == "7FFFFFFFFFFFFFFF");
}
