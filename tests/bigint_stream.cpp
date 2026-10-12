#include "cp-algo/math/bigint.hpp"
#include <cassert>
#include <sstream>
#include <type_traits>
using namespace cp_algo::math;
int main() {
    std::ostringstream out;
    auto& result = print_bigint(out, bigint<>(0));
    static_assert(std::is_same_v<decltype(result), std::ostringstream&>);
    assert(&result == &out);
    out << ' ';
    print_bigint(out, bigint<>(12));
    out << ' ';
    print_bigint(out, bigint<>(-34));
    assert(out.str() == "0 12 -34");
}
