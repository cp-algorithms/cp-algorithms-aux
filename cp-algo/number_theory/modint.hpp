#ifndef CP_ALGO_NUMBER_THEORY_MODINT_HPP
#define CP_ALGO_NUMBER_THEORY_MODINT_HPP
#include "../math/common.hpp"
#include <iostream>
#include <cassert>
namespace cp_algo::math {
    inline constexpr auto inv2(auto x) {
        assert(x % 2);
        std::make_unsigned_t<decltype(x)> y = 1;
        while(y * x != 1) {
            y *= 2 - x * y;
        }
        return y;
    }

    template<typename modint, typename _Int>
    struct modint_base {
        using Int = _Int;
        using Uint = std::make_unsigned_t<Int>;
        static constexpr size_t bits = sizeof(Int) * 8;
        using Int2 = std::conditional_t<bits <= 32, uint64_t, __uint128_t>;
        static Int mod() {
            return modint::mod();
        }
        static Uint imod() {
            return modint::imod();
        }
        static Int2 pw128() {
            return modint::pw128();
        }
        static Uint m_reduce(Int2 ab) {
            if(mod() % 2 == 0) [[unlikely]] {
                return ab % mod();
            } else {
                Uint m = ab * imod();
                return (ab + (Int2)m * mod()) >> bits;
            }
        }
        static Uint m_transform(Uint a) {
            if(mod() % 2 == 0) [[unlikely]] {
                return a;
            } else {
                return m_reduce(a * pw128());
            }
        }
        modint_base(): r(0) {}
        modint_base(Int rr): r(rr % mod()) {
            r = std::min(r, r + mod());
            r = m_transform(r);
        }
        modint inv() const {
            return bpow(to_modint(), mod() - 2);
        }
        modint operator - () const {
            modint neg;
            neg.r = std::min(-r, 2 * mod() - r);
            return neg;
        }
        modint& operator /= (const modint &t) {
            return to_modint() *= t.inv();
        }
        modint& operator *= (const modint &t) {
            r = m_reduce((Int2)r * t.r);
            return to_modint();
        }
        modint& operator += (const modint &t) {
            r += t.r; r = std::min(r, r - 2 * mod());
            return to_modint();
        }
        modint& operator -= (const modint &t) {
            r -= t.r; r = std::min(r, r + 2 * mod());
            return to_modint();
        }
        modint operator + (const modint &t) const {return modint(to_modint()) += t;}
        modint operator - (const modint &t) const {return modint(to_modint()) -= t;}
        modint operator * (const modint &t) const {return modint(to_modint()) *= t;}
        modint operator / (const modint &t) const {return modint(to_modint()) /= t;}
        // Why <=> doesn't work?..
        auto operator == (const modint_base &t) const {return getr() == t.getr();}
        auto operator != (const modint_base &t) const {return getr() != t.getr();}
        auto operator <= (const modint_base &t) const {return getr() <= t.getr();}
        auto operator >= (const modint_base &t) const {return getr() >= t.getr();}
        auto operator < (const modint_base &t) const {return getr() < t.getr();}
        auto operator > (const modint_base &t) const {return getr() > t.getr();}
        Int rem() const {
            Uint R = getr();
            return 2 * R > (Uint)mod() ? R - mod() : R;
        }

        // Only use if you really know what you're doing!
        Uint modmod() const {return (Uint)8 * mod() * mod();};
        void add_unsafe(Uint t) {r += t;}
        void pseudonormalize() {r = std::min(r, r - modmod());}
        modint const& normalize() {
            if(r >= (Uint)mod()) {
                r %= mod();
            }
            return to_modint();
        }
        void setr(Uint rr) {r = m_transform(rr);}
        Uint getr() const {
            Uint res = m_reduce(r);
            return std::min(res, res - mod());
        }
        void setr_direct(Uint rr) {r = rr;}
        Uint getr_direct() const {return r;}
    private:
        Uint r;
        modint& to_modint() {return static_cast<modint&>(*this);}
        modint const& to_modint() const {return static_cast<modint const&>(*this);}
    };
    template<typename modint>
    concept modint_type = std::is_base_of_v<modint_base<modint, typename modint::Int>, modint>;
    template<modint_type modint>
    std::istream& operator >> (std::istream &in, modint &x) {
        typename modint::Uint r;
        auto &res = in >> r;
        x.setr(r);
        return res;
    }
    template<modint_type modint>
    std::ostream& operator << (std::ostream &out, modint const& x) {
        return out << x.getr();
    }

    template<auto m>
    struct modint: modint_base<modint<m>, decltype(m)> {
        using Base = modint_base<modint<m>, decltype(m)>;
        using Base::Base;
        static constexpr Base::Uint im = m % 2 ? inv2(-m) : 0;
        static constexpr Base::Uint r2 = (typename Base::Int2)(-1) % m + 1;
        static constexpr Base::Int mod() {return m;}
        static constexpr Base::Uint imod() {return im;}
        static constexpr Base::Int2 pw128() {return r2;}
    };

    template<typename Int = int64_t>
    struct dynamic_modint: modint_base<dynamic_modint<Int>, Int> {
        using Base = modint_base<dynamic_modint<Int>, Int>;
        using Base::Base;
        static Int mod() {return m;}
        static Base::Uint imod() {return im;}
        static Base::Int2 pw128() {return r2;}
        static void switch_mod(Int nm) {
            m = nm;
            im = m % 2 ? inv2(-m) : 0;
            r2 = (typename Base::Int2)(-1) % m + 1;
        }

        // Wrapper for temp switching
        auto static with_mod(Int tmp, auto callback) {
            struct scoped {
                Int prev = mod();
                ~scoped() {switch_mod(prev);}
            } _;
            switch_mod(tmp);
            return callback();
        }
    private:
        static Int m;
        static Base::Uint im, r1, r2;
    };
    template<typename Int>
    Int dynamic_modint<Int>::m = 1;
    template<typename Int>
    dynamic_modint<Int>::Base::Uint dynamic_modint<Int>::im = -1;
    template<typename Int>
    dynamic_modint<Int>::Base::Uint dynamic_modint<Int>::r2 = 0;
}
#endif // CP_ALGO_MATH_MODINT_HPP
