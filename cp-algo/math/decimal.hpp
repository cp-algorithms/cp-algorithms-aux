#ifndef CP_ALGO_MATH_DECIMAL_HPP
#define CP_ALGO_MATH_DECIMAL_HPP
#include "bigint.hpp"
#include <utility>

namespace cp_algo::math {
    template<base_v base = x10>
    struct decimal {
        bigint<base> value;
        int64_t scale; // value * base^scale

        decimal(int64_t v=0, int64_t s=0): value(bigint<base>(v)), scale(s) {}
        decimal(bigint<base> v, int64_t s=0): value(std::move(v)), scale(s) {}

        decimal& operator *= (const decimal &other) {
            value *= other.value;
            scale += other.scale;
            return *this;
        }
        decimal& operator += (decimal const& other) {
            if (scale < other.scale) {
                value += other.value.pad(other.scale - scale);
            } else {
                value.pad_inplace(scale - other.scale);
                value += other.value;
                scale = other.scale;
            }
            return *this;
        }
        decimal& operator -= (decimal const& other) {
            if (scale < other.scale) {
                value -= other.value.pad(other.scale - scale);
            } else {
                value.pad_inplace(scale - other.scale);
                value -= other.value;
                scale = other.scale;
            }
            return *this;
        }
        decimal operator * (const decimal &other) const {
            return decimal(*this) *= other;
        }
        decimal operator + (const decimal &other) const {
            return decimal(*this) += other;
        }
        decimal operator - (const decimal &other) const {
            return decimal(*this) -= other;
        }
        auto split() const {
            auto int_part = scale >= -ssize(value.digits) ? value.top(ssize(value.digits) + scale) : bigint<base>(0);
            auto frac_part = *this - decimal(int_part);
            return std::pair{std::move(int_part), std::move(frac_part)};
        }
        void print() {
            auto [int_part, frac_part] = split();
            if(value.negative && int_part.digits.empty() && !value.digits.empty()) {
                std::cout << '-';
            }
            print_bigint(std::cout, int_part);
            if (!frac_part.value.digits.empty()) {
                std::cout << '.';
                std::cout << std::string(bigint<base>::digit_length * (-frac_part.magnitude()), '0');
                frac_part.value.negative = false;
                print_bigint<true>(std::cout, frac_part.value);
            }
            std::cout << std::endl;
        }
        bigint<base> trunc() const {
            if (scale >= 0) {
                return value.pad(scale);
            } else if (-scale >= ssize(value.digits)) {
                return 0;
            } else {
                return value.top(ssize(value.digits) + scale);
            }
        }
        bigint<base> round() const {
            if (scale >= 0) {
                return value.pad(scale);
            } else if (-scale > ssize(value.digits)) {
                return 0;
            } else {
                auto res = value.top(ssize(value.digits) + scale);
                if (value.digits[-scale - 1] * 2 >= bigint<base>::Base) {
                    res += value.negative ? -1 : 1;
                }
                return res;
            }
        }
        decimal trunc(size_t digits) const {
            digits = std::min(digits, size(value.digits));
            return decimal(
                value.top(digits),
                scale + ssize(value.digits) - digits
            );
        }
        auto magnitude() const {
            static constexpr int64_t inf = 1e18;
            if (value.digits.empty()) return -inf;
            return ssize(value.digits) + scale;
        }
        // The value must be nonzero.
        decimal inv(int64_t precision) {
            assert(precision >= 0);
            int64_t lead = llround((double)bigint<base>::Base / (double)value.digits.back());
            if(value.negative) {lead = -lead;}
            decimal d(bigint<base>(lead), -ssize(value.digits) - scale);
            size_t cur = 2;
            decimal amend = decimal(1) - trunc(cur) * d;
            // A truncated divisor can have an exact reciprocal while its omitted
            // low limbs still matter. Include enough input limbs before stopping.
            while(cur <= size_t(precision) || -amend.magnitude() < precision) {
                d += d * amend;
                cur = 2 * (1 - amend.magnitude());
                d = d.trunc(cur);
                amend = decimal(1) - trunc(cur) * d;
            }
            return d;
        }
    };

    template<base_v base>
    auto divmod_fast(bigint<base> const& a, int64_t b) {
        // Optimized divmod for small divisors that fit in int64_t
        if (b == 0) {
            assert(false && "Division by zero");
        }
        bool neg_a = a.negative;
        bool neg_b = b < 0;
        uint64_t divisor = neg_b ? -uint64_t(b) : uint64_t(b);
        
        bigint<base> quotient;
        uint64_t remainder = 0;

        auto n = ssize(a.digits);
        quotient.digits.resize(n);
        for (auto i = n - 1; i >= 0; i--) {
            __uint128_t val = (__uint128_t)remainder * bigint<base>::Base + a.digits[i];
            uint64_t q = uint64_t(val / divisor);
            remainder = uint64_t(val % divisor);
            quotient.digits[i] = q;
        }
        quotient.negative = (neg_a ^ neg_b);
        quotient.normalize();
        
        bigint<base> rem{int64_t(remainder)};
        rem.negative = neg_a && remainder != 0;
        
        return std::pair{std::move(quotient), std::move(rem)};
    }

    // Truncate the quotient toward zero; the divisor must be nonzero.
    template<base_v base>
    std::pair<bigint<base>, bigint<base>> divmod(bigint<base> const& a, bigint<base> const& b) {
        if(a.negative || b.negative) {
            auto magnitude_a = a, magnitude_b = b;
            magnitude_a.negative = magnitude_b.negative = false;
            auto [quotient, remainder] = divmod(magnitude_a, magnitude_b);
            quotient.negative = a.negative != b.negative;
            remainder.negative = a.negative;
            quotient.normalize();
            remainder.normalize();
            return std::pair{std::move(quotient), std::move(remainder)};
        }
        if (a < b) {
            return std::pair{bigint<base>(0), a};
        }
        // Use fast path if b fits in int64_t
        if (size(b.digits) == 1) {
            int64_t b_val = b.negative ? -int64_t(b.digits[0]) : int64_t(b.digits[0]);
            return divmod_fast(a, b_val);
        }
        // General case using decimal arithmetic
        auto A = decimal<base>(a);
        auto B = decimal<base>(b);
        auto d = (A * B.inv(A.magnitude() - B.magnitude() + 1)).round();
        auto r = a - d * b;
        if (r >= b) {
            d += 1;
            r -= b;
        }
        if (r < bigint<base>(0)) {
            d -= 1;
            r += b;
        }
        return std::pair{std::move(d), std::move(r)};
    }
}

#endif // CP_ALGO_MATH_DECIMAL_HPP
