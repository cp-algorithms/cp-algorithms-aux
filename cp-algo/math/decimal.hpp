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
        decimal(bigint<base> v, int64_t s=0): value(v), scale(s) {}

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
            return std::pair{int_part, frac_part};
        }
        void print() {
            auto [int_part, frac_part] = split();
            print_bigint(std::cout, int_part);
            if (frac_part.value != bigint<base>(0)) {
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
                    res += 1;
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
        decimal inv(int64_t precision) {
            assert(precision >= 0);
            int64_t lead = llround((double)bigint<base>::Base / (double)value.digits.back());
            decimal d(bigint<base>(lead), -ssize(value.digits));
            size_t cur = 2;
            decimal amend = decimal(1) - trunc(cur) * d;
            while(-amend.magnitude() < precision) {
                d += d * amend;
                cur = 2 * (1 - amend.magnitude());
                d = d.trunc(cur);
                amend = decimal(1) - trunc(cur) * d;
            }
            return d;
        }
    };

    template<base_v base>
    auto divmod(bigint<base> const& a, bigint<base> const& b) {
        if (a < b) {
            return std::pair{bigint<base>(0), a};
        }
        auto A = decimal<base>(a);
        auto B = decimal<base>(b);
        auto d = (A * B.inv(A.magnitude())).trunc();
        auto r = a - d * b;
        if (r >= b) {
            d += 1;
            r -= b;
        }
        return std::pair{d, r};
    }
}

#endif // CP_ALGO_MATH_DECIMAL_HPP
