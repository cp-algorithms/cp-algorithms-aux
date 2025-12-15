#ifndef CP_ALGO_MATH_BIGINT_HPP
#define CP_ALGO_MATH_BIGINT_HPP
#include "../util/big_alloc.hpp"
#include "../math/fft_simple.hpp"
#include <bits/stdc++.h>

namespace cp_algo::math {
    enum base_v {
        x10 = uint64_t(1e16),
        x16 = uint64_t(1ull << 60)
    };
    template<base_v base = x10>
    struct bigint {
        static constexpr uint64_t Base = uint64_t(base);
        static constexpr uint16_t digit_length = base == x10 ? 16 : 15;
        static constexpr uint16_t sub_base = base == x10 ? 10 : 16;
        static constexpr uint32_t meta_base = base == x10 ? uint32_t(1e4) : uint32_t(1 << 15);
        big_basic_string<uint64_t> digits;
        bool negative;

        auto operator <=> (bigint const& other) const {
            // Handle zero cases
            if (digits.empty() && other.digits.empty()) {
                return std::strong_ordering::equal;
            }
            if (digits.empty()) {
                return other.negative ? std::strong_ordering::greater : std::strong_ordering::less;
            }
            if (other.digits.empty()) {
                return negative ? std::strong_ordering::less : std::strong_ordering::greater;
            }
            
            // Handle sign differences
            if (negative != other.negative) {
                return negative ? std::strong_ordering::less : std::strong_ordering::greater;
            }
            
            // Both have the same sign - compare magnitudes
            if (digits.size() != other.digits.size()) {
                auto size_cmp = digits.size() <=> other.digits.size();
                // If both negative, reverse the comparison
                return negative ? 0 <=> size_cmp : size_cmp;
            }
            
            // Same size, compare digits from most significant to least
            for (auto i = ssize(digits) - 1; i >= 0; i--) {
                auto digit_cmp = digits[i] <=> other.digits[i];
                if (digit_cmp != std::strong_ordering::equal) {
                    return negative ? 0 <=> digit_cmp : digit_cmp;
                }
            }
            
            return std::strong_ordering::equal;
        }

        bigint() {}

        bigint(big_basic_string<uint64_t> d, bool neg): digits(std::move(d)), negative(neg) {
            normalize();
        }

        bigint& pad_inplace(size_t to_add) {
            digits.insert(0, to_add, 0);
            return normalize();
        }
        bigint& drop_inplace(size_t to_drop) {
            digits.erase(0, std::min(to_drop, size(digits)));
            return normalize();
        }
        bigint& take_inplace(size_t to_keep) {
            digits.erase(std::min(to_keep, size(digits)), std::string::npos);
            return normalize();
        }
        bigint& top_inplace(size_t to_keep) {
            if (to_keep >= size(digits)) {
                return pad_inplace(to_keep - size(digits));
            } else {
                return drop_inplace(size(digits) - to_keep);
            }
        }
        bigint pad(size_t to_add) const {
            return bigint{big_basic_string<uint64_t>(to_add, 0) + digits, negative}.normalize();
        }
        bigint drop(size_t to_drop) const {
            return bigint{digits.substr(std::min(to_drop, size(digits))), negative}.normalize();
        }
        bigint take(size_t to_keep) const {
            return bigint{digits.substr(0, std::min(to_keep, size(digits))), negative}.normalize();
        }
        bigint top(size_t to_keep) const {
            if (to_keep >= size(digits)) {
                return pad(to_keep - size(digits));
            } else {
                return drop(size(digits) - to_keep);
            }
        }

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
                d_ptr[i] -= o_ptr[i] + carry;
                carry = d_ptr[i] >= base;
                d_ptr[i] += carry ? uint64_t(base) : 0;
            }
            if (carry) {
                N = size(digits);   
                for (; i < N && d_ptr[i] == 0; i++) {
                    d_ptr[i] = base - 1;
                }
                if (i < N) {
                    d_ptr[i]--;
                } else {
                    // Two's complement: flip all digits then add 1
                    for (i = 0; i < N; i++) {
                        d_ptr[i] = base - d_ptr[i] - 1;
                    }
                    bool carry = true;
                    for (i = 0; i < N && carry; i++) {
                        d_ptr[i]++;
                        carry = d_ptr[i] >= base;
                        d_ptr[i] -= carry * base;
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
                d_ptr[i] += o_ptr[i] + carry;
                carry = d_ptr[i] >= base;
                d_ptr[i] -= carry ? uint64_t(base) : 0;
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
        bigint(int64_t x) {
            negative = x < 0;
            x = negative ? -x : x;
            digits = x ? big_basic_string<uint64_t>{uint64_t(x)} : big_basic_string<uint64_t>{};
        }
        bigint(std::span<char> s): negative(false) {
            if (size(s) < digit_length) {
                int64_t val = 0;
                std::from_chars(s.data(), s.data() + size(s), val, sub_base);
                *this = bigint(val);
                return;
            }
            if (!empty(s) && s[0] == '-') {
                negative = true;
                s = s.subspan(1);
            }
            size_t len = size(s);
            assert(len > 0);
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
        void to_metabase() {
            auto N = ssize(digits);
            digits.resize(4 * N);
            for (auto i = N - 1; i >= 0; i--) {
                uint64_t val = digits[i];
                digits[4 * i] = val % meta_base;
                val /= meta_base;
                digits[4 * i + 1] = val % meta_base;
                val /= meta_base;
                digits[4 * i + 2] = val % meta_base;
                val /= meta_base;
                digits[4 * i + 3] = val;
            }
        }
        void from_metabase() {
            auto N = (ssize(digits) + 3) / 4;
            digits.resize(4 * N);
            uint64_t carry = 0;
            for (int i = 0; i < N; i++) {
                __uint128_t val = digits[4 * i + 3];
                val = val * meta_base + digits[4 * i + 2];
                val = val * meta_base + digits[4 * i + 1];
                val = val * meta_base + digits[4 * i];
                val += carry;
                digits[i] = uint64_t(val % base);
                carry = uint64_t(val / base);
            }
            digits.resize(N);
            while (carry) {
                digits.push_back(carry % base);
                carry /= base;
            }
        }
        bigint& operator *= (int64_t other) {
            if (other < 0) {
                negative ^= 1;
                other = -other;
            }
            if (other == 0) {
                return *this = bigint(0);
            } else if (other == 1) {
                return *this;
            }
            uint64_t carry = 0;
            for (auto &d: digits) {
                __uint128_t val = __uint128_t(d) * other + carry;
                d = uint64_t(val % base);
                carry = uint64_t(val / base);
            }
            if (carry) {
                digits.push_back(carry % base);
                carry /= base;
            }
            return *this;
        }
        bigint operator * (int64_t other) const {
            return bigint(*this) *= other;
        }
        friend bigint operator * (int64_t lhs, const bigint& rhs) {
            return bigint(rhs) *= lhs;
        }
        bigint& mul_inplace(auto &&other) {
            negative ^= other.negative;
            auto n = size(digits), m = size(other.digits);
            if (n < m) {
                std::swap(n, m);
                std::swap(digits, other.digits);
            }
            if (m <= 1) {
                return *this *= int64_t(m == 0 ? 0 : other.digits[0]);
            }
            to_metabase();
            other.to_metabase();
            fft::conv_simple(digits, other.digits);
            from_metabase();
            return normalize();
        }
        bigint& operator *= (bigint const& other) {
            return mul_inplace(bigint(other));
        }
        bigint operator * (const bigint& other) const {
            return bigint(*this).mul_inplace(bigint(other));
        }
    };

    template<base_v base>
    decltype(std::cin)& operator >> (decltype(std::cin) &in, cp_algo::math::bigint<base> &x) {
        std::string s;
        in >> s;
        x = {s};
        return in;
    }

    template<base_v base, bool fill = true>
    auto& print_digit(auto &out, uint64_t d) {
        char buf[16];
        auto [ptr, ec] = std::to_chars(buf, buf + sizeof(buf), d, bigint<base>::sub_base);
        if constexpr (base == x16) {
            std::ranges::transform(buf, buf, toupper);
        }
        auto len = ptr - buf;
        if constexpr (fill) {
            out << std::string(bigint<base>::digit_length - len, '0');
        }
        return out << std::string_view(buf, len);
    }

    template<bool fill_all = false, base_v base>
    auto& print_bigint(auto &out, cp_algo::math::bigint<base> const& x) {
        if (x.negative) {
            out << '-';
        }
        if (empty(x.digits)) {
            return print_digit<base, fill_all>(out, 0);
        }
        print_digit<base, fill_all>(out, x.digits.back());
        for (auto d: x.digits | std::views::reverse | std::views::drop(1)) {
            print_digit<base, true>(out, d);
        }
        return out;
    }

    template<base_v base>
    decltype(std::cout)& operator << (decltype(std::cout) &out, cp_algo::math::bigint<base> const& x) {
        return print_bigint(out, x);
    }
}

#endif // CP_ALGO_MATH_BIGINT_HPP
