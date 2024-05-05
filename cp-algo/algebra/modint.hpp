#ifndef CP_ALGO_ALGEBRA_MODINT_HPP
#define CP_ALGO_ALGEBRA_MODINT_HPP
#include "common.hpp"
#include <iostream>
namespace cp_algo::algebra {
    template<typename modint>
    struct modint_base {
        static int mod() {
            return modint::mod();
        }
        modint_base(): r(0) {}
        modint_base(int64_t rr): r(rr % mod()) {
            r = std::min(r, r + mod());
        }
        modint inv() const {
            return bpow(to_modint(), mod() - 2);
        }
        modint operator - () const {return std::min(-r, mod() - r);}
        modint& operator /= (const modint &t) {
            return to_modint() *= t.inv();
        }
        modint& operator *= (const modint &t) {
            r *= t.r; if(mod()) {r %= mod();}
            return to_modint();
        }
        modint& operator += (const modint &t) {
            r += t.r; r = std::min(r, r - mod());
            return to_modint();
        }
        modint& operator -= (const modint &t) {
            r -= t.r; r = std::min(r, r + mod());
            return to_modint();
        }
        modint operator + (const modint &t) const {return modint(to_modint()) += t;}
        modint operator - (const modint &t) const {return modint(to_modint()) -= t;}
        modint operator * (const modint &t) const {return modint(to_modint()) *= t;}
        modint operator / (const modint &t) const {return modint(to_modint()) /= t;}
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
            return to_modint();
        }
        uint64_t& setr() {return r;}
        uint64_t getr() const {return r;}
    private:
        uint64_t r;
        modint& to_modint() {return static_cast<modint&>(*this);}
        modint const& to_modint() const {return static_cast<modint const&>(*this);}
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
        static constexpr int mod() {return m;}
        using Base = modint_base<modint<m>>;
        using Base::Base;
    };

    struct dynamic_modint: modint_base<dynamic_modint> {
        static int mod() {return m;}
        static void switch_mod(int nm) {m = nm;}
        using Base = modint_base<dynamic_modint>;
        using Base::Base;
    private:
        static int m;
    };
    int dynamic_modint::m = 0;
}
#endif // CP_ALGO_ALGEBRA_MODINT_HPP
