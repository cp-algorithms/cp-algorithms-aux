#include "cp-algo/math/bigint.hpp"
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

using namespace cp_algo::math;

void check(size_t n, size_t m) {
    assert(n >= m);
    std::string sa(n, 'F'), sb(m, 'F');
    bigint<x16> a(sa), b(sb);
    auto product = a * b;
    std::ostringstream out;
    out << product;
    // (16^n - 1)(16^m - 1), including the long borrow/carry chains.
    std::string expected = std::string(m - 1, 'F') + 'E'
                         + std::string(n - m, 'F') + std::string(m - 1, '0') + '1';
    assert(out.str() == expected);
}

void check_sparse(size_t n, size_t m) {
    std::string sa(n + 1, '0'), sb(m + 1, '0');
    sa.front() = sa.back() = sb.front() = sb.back() = '1';
    bigint<x16> a(sa), b(sb);
    std::ostringstream out;
    out << a * b;
    std::string expected(n + m + 1, '0');
    for(size_t exponent: {size_t(0), n, m, n + m}) {
        expected[n + m - exponent]++;
    }
    assert(out.str() == expected);
}

void check_balanced_boundary() {
    constexpr size_t n = 262144;
    constexpr uint64_t digit = 16384;
    uint64_t limb = 0;
    for(int j = 0; j < 4; j++) limb |= digit << (15 * j);
    bigint<x16> a(cp_algo::big_basic_string<uint64_t>(n / 4, limb), false), b = a;
    auto product = a * b;
    cp_algo::big_basic_string<uint64_t> expected(n / 2, 0);
    uint64_t carry = 0;
    for(size_t k = 0; k < 2 * n; k++) {
        uint64_t count = k < n ? k + 1 : 2 * n - 1 - k;
        uint64_t value = count * digit * digit + carry;
        expected[k / 4] |= (value & 32767) << (15 * (k % 4));
        carry = value >> 15;
    }
    assert(!carry && product.digits == expected && !product.negative);
    a.negate();
    auto negative = a * b;
    assert(negative.negative && negative.digits == expected);
    b.negate();
    auto positive = a * b;
    assert(!positive.negative && positive.digits == expected);
    bigint<x16> zero(0);
    zero.negate();
    auto result = a * zero;
    assert(result.digits.empty() && !result.negative);
}

int main() {
    // Both sides of the schoolbook/FFT cutoff and partial 60-bit limbs.
    check(480, 480);
    check(495, 495);
    check(511, 496);
    // Uncentered 15-bit digits round individual convolution coefficients
    // incorrectly on these supported sizes, corrupting the following carries.
    check(1504000, 1504000);
    check(1565000, 1565000);
    check(1600000, 1600000);
    check(1600000, 800000);
    check(1600000, 497);
    // Centering introduces large negative FFT digits even for sparse inputs.
    // Their exact correction must recover long runs of zero coefficients.
    check_sparse(1599999, 1599999);
    check_sparse(1599999, 799999);
    // Half-base digits require negative FFT digits and signed floor carries.
    // This length also sits at a transform boundary: no extra digit is needed.
    check_balanced_boundary();
    std::cout << "Hex multiplication FFT precision and carry checks passed\n";
}
