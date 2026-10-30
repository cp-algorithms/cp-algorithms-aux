// Licensed under GPL-3.0-only with the online-judge exception.
// See cp-algo/LICENSE and cp-algo/ONLINE-JUDGE-EXCEPTION.
#ifndef CP_ALGO_NUMBER_THEORY_MODINT_HPP
#define CP_ALGO_NUMBER_THEORY_MODINT_HPP
#include "../math/common.hpp"
#include <iostream>
#include <cassert>
#include <concepts>
#include <type_traits>
namespace cp_algo::math {

    // M > 0 fixes the modulus; M == 0 uses a thread-local runtime modulus.
    template<auto M>
    struct modint {
        using Int = decltype(M);
        using UInt = std::make_unsigned_t<Int>;
        static constexpr size_t bits = sizeof(Int) * 8;
        using Int2 = std::conditional_t<bits <= 32, int64_t, __int128_t>;
        using UInt2 = std::conditional_t<bits <= 32, uint64_t, __uint128_t>;
        constexpr static Int mod() {
            if constexpr(M == 0) {return m;}
            else {return M;}
        }
        constexpr static UInt remod() {
            if constexpr(M == 0) {return rm;}
            else {return M;}
        }
        constexpr static UInt2 modmod() {
            return UInt2(mod()) * mod();
        }
        constexpr modint() = default;
        constexpr modint(Int2 rr) {
            setr(UInt((rr + modmod()) % mod()));
        }
        constexpr modint inv() const {
            return bpow(*this, mod() - 2);
        }
        modint operator-() const {
            modint neg;
            neg.r = std::min(-r, remod() - r);
            return neg;
        }
        modint& operator/=(const modint& t) {
            return *this *= t.inv();
        }
        modint& operator*=(const modint& t) {
            if constexpr(M == 0) {r = m_reduce(UInt2(r) * t.r);}
            else {r = UInt(UInt2(r) * t.r % mod());}
            return *this;
        }
        modint& operator+=(const modint& t) {
            r += t.r;
            r = std::min(r, r - remod());
            return *this;
        }
        modint& operator-=(const modint& t) {
            r -= t.r;
            r = std::min(r, r + remod());
            return *this;
        }
        modint operator+(const modint& t) const { return modint(*this) += t; }
        modint operator-(const modint& t) const { return modint(*this) -= t; }
        modint operator*(const modint& t) const { return modint(*this) *= t; }
        modint operator/(const modint& t) const { return modint(*this) /= t; }
        // Why <=> doesn't work?..
        auto operator==(const modint& t) const { return getr() == t.getr(); }
        auto operator!=(const modint& t) const { return getr() != t.getr(); }
        auto operator<=(const modint& t) const { return getr() <= t.getr(); }
        auto operator>=(const modint& t) const { return getr() >= t.getr(); }
        auto operator<(const modint& t) const { return getr() < t.getr(); }
        auto operator>(const modint& t) const { return getr() > t.getr(); }
        Int rem() const {
            UInt R = getr();
            return R - (R > (UInt)mod() / 2) * mod();
        }
        constexpr void setr(UInt rr) {
            if constexpr(M == 0) {r = m_transform(rr);}
            else {r = rr;}
        }
        constexpr UInt getr() const {
            if constexpr(M == 0) {
                UInt res = m_reduce(r);
                return std::min(res, res - mod());
            } else {return r;}
        }
        // Only use these if you really know what you're doing!
        static uint64_t modmod8() { return uint64_t(8 * modmod()); }
        void add_unsafe(UInt t) { r += t; }
        void pseudonormalize() { r = std::min(r, r - modmod8()); }
        modint const& normalize() {
            if(r >= (UInt)mod()) {
                r %= mod();
            }
            return *this;
        }
        void setr_direct(UInt rr) { r = rr; }
        UInt getr_direct() const { return r; }

        // Odd moduli up to a quarter of the unsigned word keep Montgomery residues lazily in
        // [0, 2 mod): remod() = 2 mod, and both 4 mod and ab + q * mod fit. Any other modulus keeps
        // fully reduced residues, remod() = mod, which is what the sums of modint need then:
        // an even one without the Montgomery form, a wide odd one with a reduction that subtracts
        // high words instead of adding double words that would overflow.
        // Out of line, so that the hot path stays as small as it was.
        [[gnu::noinline, gnu::cold]] static UInt m_reduce_reduced(UInt2 ab) requires(M == 0) {
            if(mod() % 2 == 0) {
                return UInt(ab % mod());
            }
            // q * mod has the low word of ab, so the difference of the high words is exact.
            UInt q = -(UInt(ab) * inverse);
            auto high = UInt(ab >> bits);
            auto low = UInt(UInt2(q) * UInt(mod()) >> bits);
            return high >= low ? high - low : high - low + mod();
        }
        static UInt m_reduce(UInt2 ab) requires(M == 0) {
            if(imod() == 0) [[unlikely]] {
                return m_reduce_reduced(ab);
            }
            UInt2 m = UInt(ab) * imod();
            return UInt((ab + m * mod()) >> bits);
        }
        static UInt m_transform(UInt a) requires(M == 0) {
            if(mod() % 2 == 0) [[unlikely]] {
                return a;
            } else {
                return m_reduce(a * pw128());
            }
        }
        static UInt imod() requires(M == 0) { return im; }
        static UInt2 pw128() requires(M == 0) { return r2; }
        static void switch_mod(Int nm) requires(M == 0) {
            m = nm;
            bool lazy = m % 2 && UInt(m) <= UInt(-1) / 4;
            rm = UInt(m) * (lazy ? 2 : 1);
            inverse = m % 2 ? inv2(-m) : 0;
            im = lazy ? inverse : 0;
            r2 = static_cast<UInt>(static_cast<UInt2>(-1) % m + 1);
        }

        // Wrapper for temp switching
        auto static with_mod(Int tmp, auto callback) requires(M == 0) {
            struct scoped {
                Int prev = mod();
                ~scoped() { switch_mod(prev); }
            } _;
            switch_mod(tmp);
            return callback();
        }

    private:
        UInt r;
        static inline thread_local Int m = 1;
        // im: -1 / mod modulo 2^bits for lazy residues and 0 for reduced ones; inverse: the same
        // for every odd mod; rm: the value of remod().
        static inline thread_local UInt im = -1, r2 = 0, inverse = -1, rm = 2;
    };

    template<typename Int = int>
    using dynamic_modint = modint<Int(0)>;

    namespace impl {
        template<typename T> constexpr bool is_modint = false;
        template<auto M> constexpr bool is_modint<modint<M>> = true;
    }
    template<typename T>
    concept modint_type = impl::is_modint<T>;
    // Choose the input integer type explicitly. Reduce before narrowing to storage.
    // Failed extraction leaves x unchanged, including for signed minima.
    template<std::integral Int, auto M, class Stream> requires (!std::same_as<Int, bool>)
    void read(modint<M>& x, Stream& in, std::type_identity<Int>) {
        using base = modint<M>;
        using UInt = std::common_type_t<std::make_unsigned_t<Int>, typename base::UInt>;
        Int value;
        if(in >> value) {
            // Convert before negating so the signed minimum is supported too.
            UInt magnitude = value < 0 ? -UInt(value) : UInt(value);
            bool reduce = true;
            if constexpr(M == 0 && sizeof(UInt) == sizeof(typename base::UInt)) {
                // Montgomery conversion already reduces any word for odd moduli.
                reduce = base::mod() % 2 == 0;
            }
            if(reduce && magnitude >= UInt(base::mod())) {magnitude %= UInt(base::mod());}
            x.setr(typename base::UInt(magnitude));
            if(value < 0) {x = -x;}
        }
    }
    // Input must already satisfy 0 <= value < mod(); no reduction or range check.
    // Failed extraction leaves x unchanged.
    template<auto M, class Stream>
    void read_reduced(modint<M>& x, Stream& in, std::type_identity<void>) {
        typename modint<M>::UInt value;
        if(in >> value) {x.setr(value);}
    }
    template<auto M, class Stream>
    void write(modint<M> const& x, Stream& out, std::type_identity<void>) {
        out << x.getr();
    }
    template<modint_type modint>
    decltype(std::cout)& operator<<(decltype(std::cout)& out, modint const& x) {
        return out << x.getr();
    }

}
#endif // CP_ALGO_MATH_MODINT_HPP
