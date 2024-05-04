#ifndef CP_ALGO_ALGEBRA_MODINT_HPP
#define CP_ALGO_ALGEBRA_MODINT_HPP
#include "../random/rng.hpp"
#include "affine.hpp"
#include "common.hpp"
#include <algorithm>
#include <iostream>
#include <optional>
namespace cp_algo::algebra {
    template<typename modint>
    struct modint_base {
        // Would make it virtual, but it affects avx...
        int mod() const {
            if constexpr(modint::static_mod) {
                return modint::static_mod;
            } else {
                return static_cast<modint*>(this)->mod();
            }
        }
        modint_base(): r(0) {}
        modint_base(int64_t rr): r(rr % mod()) {
            r = std::min(r, r + mod());
        }
        modint inv() const {
            return bpow(static_cast<modint const&>(*this), mod() - 2);
        }
        modint operator - () const {return std::min(-r, mod() - r);}
        modint& operator /= (const modint &t) {return *this *= t.inv();}
        modint& operator *= (const modint &t) {
            r *= t.r; r %= mod();
            return static_cast<modint&>(*this);
        }
        modint& operator += (const modint &t) {
            r += t.r; r = std::min(r, r - mod());
            return static_cast<modint&>(*this);
        }
        modint& operator -= (const modint &t) {
            r -= t.r; r = std::min(r, r + mod());
            return static_cast<modint&>(*this);
        }
        modint operator + (const modint &t) const {
            return modint(static_cast<modint const&>(*this)) += t;
        }
        modint operator - (const modint &t) const {
            return modint(static_cast<modint const&>(*this)) -= t;
        }
        modint operator * (const modint &t) const {
            return modint(static_cast<modint const&>(*this)) *= t;
        }
        modint operator / (const modint &t) const {
            return modint(static_cast<modint const&>(*this)) /= t;
        }
        auto operator <=> (const modint_base &t) const = default;
        explicit operator int() const {return r;}
        int64_t rem() const {return 2 * r > (uint64_t)mod() ? r - mod() : r;}

        // Only use if you really know what you're doing!
        uint64_t modmod() const {return 8LL * mod() * mod();};
        void add_unsafe(uint64_t t) {r += t;}
        void pseudonormalize() {r = std::min(r, r - modmod());}
        modint const& normalize() {
            if(r >= (uint64_t)mod()) {
                r %= mod();
            }
            return static_cast<modint&>(*this);
        }
        uint64_t& setr() {return r;}
        uint64_t getr() const {return r;}
    private:
        uint64_t r;
    };
    template<typename modint>
    std::istream& operator >> (std::istream &in, modint_base<modint> &x) {

        return in >> x.setr();
    }
    template<typename modint>
    std::ostream& operator << (std::ostream &out, modint_base<modint> const& x) {
        return out << x.getr();
    }

    template<int m>
    struct modint: modint_base<modint<m>> {
        constexpr static int static_mod = m;
        using Base = modint_base<modint<m>>;
        using Base::Base;
        constexpr static uint64_t static_modmod = 8LL*m*m;
    };

    struct dynamic_modint: modint_base<dynamic_modint> {
        constexpr static int static_mod = 0;
        using Base = modint_base<dynamic_modint>;
        int mod() const {return m;}
        dynamic_modint(dynamic_modint const& t): Base(t), m(t.m) {}
        dynamic_modint(int m, int64_t r): m(m) {
            setr() = r % m;
            setr() = std::min(setr(), setr() + mod());
        }
        static auto GF(int mod) {
            return [mod](int64_t r) {
                return dynamic_modint(mod, r);
            };
        }
    private:
        int m;
    };
}
#endif // CP_ALGO_ALGEBRA_MODINT_HPP
