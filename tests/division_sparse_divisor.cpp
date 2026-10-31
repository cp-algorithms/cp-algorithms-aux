#include "cp-algo/math/decimal.hpp"
#include <cassert>
#include <sstream>
#include <string>

using namespace cp_algo::math;

template<base_v base>
void check() {
    constexpr int width = bigint<base>::digit_length;
    constexpr char high = base == x10 ? '5' : '8';
    constexpr char max = base == x10 ? '9' : 'F';
    for(int m : {2 * width - 1, 2 * width, 2 * width + 1, 4 * width + 1}) {
        // B = (radix / 2) * radix^(m-1) + 1, Q = radix^(3m) - 1.
        // Write A = B * Q directly, independently of bigint multiplication.
        std::string a = high + std::string(m - 1, '0') + std::string(2 * m, max)
                      + char(high - 1) + std::string(m - 1, max);
        std::string b = high + std::string(m - 2, '0') + '1';
        for(bool negative : {false, true}) {
            auto signed_a = (negative ? "-" : "") + a;
            auto [q, r] = divmod(bigint<base>(signed_a), bigint<base>(b));
            std::ostringstream out;
            out << q;
            assert(out.str() == (negative ? "-" : "") + std::string(3 * m, max));
            assert(r.digits.empty() && !r.negative);
        }
    }
}

int main() {
    check<x10>();
    check<x16>();
}
