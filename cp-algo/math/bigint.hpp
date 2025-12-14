#ifndef CP_ALGO_MATH_BIGINT_HPP
#define CP_ALGO_MATH_BIGINT_HPP
#include "../util/big_alloc.hpp"
#include <bits/stdc++.h>

namespace cp_algo::math {
    enum base_v {
        x10 = uint64_t(1e18),
        x16 = uint64_t(0)
    };
    template<base_v base = x10>
    struct bigint {
        big_vector<uint64_t> digits;
        bool negative;

        bigint() {}

        bigint& normalize() {
            while (!empty(digits) && digits.back() == 0) {
                digits.pop_back();
            }
            if (digits.empty()) {
                negative = false;
            }
            return *this;
        }
        bigint& negate() {
            negative ^= 1;
            return *this;
        }
        bigint operator -() {
            return bigint(*this).negate();
        }
        bigint& operator -= (const bigint& other) {
            if (negative != other.negative) {
                return (negate() += other).negate().normalize();
            }
            digits.resize(std::max(size(digits), size(other.digits)));
            bool carry = false;
            auto d_ptr = std::assume_aligned<32>(digits.data());
            auto o_ptr = std::assume_aligned<32>(other.digits.data());
            size_t N = size(other.digits);
            size_t i = 0;
            for (; i < N; i++) {
                if constexpr (base == x10) {
                    d_ptr[i] -= o_ptr[i] + carry;
                    carry = d_ptr[i] >= base;
                    d_ptr[i] += carry ? uint64_t(base) : 0;
                } else if constexpr (base == x16) {
                    auto sub = o_ptr[i] + carry;
                    auto new_carry = sub ? (d_ptr[i] < sub) : carry;
                    d_ptr[i] -= sub;
                    carry = new_carry;
                } else {
                    static_assert(base == x10 || base == x16, "Unsupported base");
                }
            }
            if (carry) {
                N = size(digits);   
                for (; i < N && d_ptr[i] == 0; i++) {
                    d_ptr[i] = base - 1;
                }
                if (i < N) {
                    d_ptr[i]--;
                } else {
                    d_ptr[0]--;
                    for (i = 0; i < N; i++) {
                        d_ptr[i] = base - d_ptr[i] - 1;
                    }
                    negate();
                }
            }
            return normalize();
        }
        bigint& operator += (const bigint& other) {
            if (negative != other.negative) {
                return (negate() -= other).negate().normalize();
            }
            digits.resize(std::max(size(digits), size(other.digits)));
            bool carry = false;
            auto d_ptr = std::assume_aligned<32>(digits.data());
            auto o_ptr = std::assume_aligned<32>(other.digits.data());
            size_t N = size(other.digits);
            size_t i = 0;
            for (; i < N; i++) {
                if constexpr (base == x10) {
                    d_ptr[i] += o_ptr[i] + carry;
                    carry = d_ptr[i] >= base;
                    d_ptr[i] -= carry ? uint64_t(base) : 0;
                } else if constexpr (base == x16) {
                    auto add = o_ptr[i] + carry;
                    auto new_carry = add ? (d_ptr[i] >= -add) : carry;
                    d_ptr[i] += add;
                    carry = new_carry;
                } else {
                    static_assert(base == x10 || base == x16, "Unsupported base");
                }
            }
            if (carry) {
                N = size(digits);
                for (; i < N && d_ptr[i] == uint64_t(base) - 1; i++) {
                    d_ptr[i] = 0;
                }
                if (i < N) {
                    d_ptr[i]++;
                } else {
                    digits.push_back(1);
                }
            }
            return *this;
        }

        bigint(std::span<char> s): negative(false) {
            if (!empty(s) && s[0] == '-') {
                negative = true;
                s = s.subspan(1);
            }
            size_t len = size(s);
            assert(len > 0);
            constexpr auto digit_length = base == x10 ? 18 : base == x16 ? 16 : 0;
            constexpr auto sub_base = base == x10 ? 10 : base == x16 ? 16 : 0;
            size_t num_digits = (len + digit_length - 1) / digit_length;
            digits.resize(num_digits);
            size_t i = len;
            for (size_t j = 0; j < num_digits - 1; j++) {
                std::from_chars(s.data() + i - digit_length, s.data() + i, digits[j], sub_base);
                i -= digit_length;
            }
            std::from_chars(s.data(), s.data() + i, digits.back(), sub_base);
            normalize();
        }

        bigint operator + (const bigint& other) const {
            return bigint(*this) += other;
        }
        bigint operator - (const bigint& other) const {
            return bigint(*this) -= other;
        }
    };

    template<base_v base>
    decltype(std::cin)& operator >> (decltype(std::cin) &in, cp_algo::math::bigint<base> &x) {
        std::string s;
        in >> s;
        x = {s};
        return in;
    }

    template<base_v base>
    decltype(std::cout)& operator << (decltype(std::cout) &out, cp_algo::math::bigint<base> const& x) {
        if (x.negative) {
            out << '-';
        }
        if (empty(x.digits)) {
            return out << '0';
        }
        constexpr auto digit_length = base == x10 ? 18 : base == x16 ? 16 : 0;
        constexpr auto sub_base = base == x10 ? 10 : base == x16 ? 16 : 0;
        char buf[20];
        auto [ptr, ec] = std::to_chars(buf, buf + sizeof(buf), x.digits.back(), sub_base);
        if constexpr (base == x16) {
            for (char* p = buf; p != ptr; ++p) {
                *p = (char)std::toupper(*p);
            }
        }
        buf[ptr - buf] = 0;
        out << buf;
        for (auto d: x.digits | std::views::reverse | std::views::drop(1)) {
            auto [ptr, ec] = std::to_chars(buf, buf + sizeof(buf), d, sub_base);
            if constexpr (base == x16) {
                for (char* p = buf; p != ptr; ++p) {
                    *p = (char)std::toupper(*p);
                }
            }
            auto len = ptr - buf;
            out << std::string(digit_length - len, '0');
            buf[len] = 0;
            out << buf;
        }
        return out;
    }
}

#endif // CP_ALGO_MATH_BIGINT_HPP
