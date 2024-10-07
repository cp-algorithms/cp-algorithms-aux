#ifndef CP_ALGO_NUMBER_THEORY_MODINT_HPP
#define CP_ALGO_NUMBER_THEORY_MODINT_HPP
#include "../math/common.hpp"
#include <iostream>
#include <cassert>
namespace cp_algo::math {
    inline constexpr uint64_t inv64(uint64_t x) {
        assert(x % 2);
        uint64_t y = 1;
        while(y * x != 1) {
            y *= 2 - x * y;
        }
        return y;
    }

    template<typename modint>
    struct modint_base {
        static int64_t mod() {
            return modint::mod();
        }
        static uint64_t imod() {
            return modint::imod();
        }
        static __uint128_t pw128() {
            return modint::pw128();
        }
        static uint64_t m_reduce(__uint128_t ab) {
            if(mod() % 2 == 0) [[unlikely]] {
                return ab % mod();
            } else {
                uint64_t m = ab * imod();
                return (ab + __uint128_t(m) * mod()) >> 64;
            }
        }
        static uint64_t m_transform(uint64_t a) {
            if(mod() % 2 == 0) [[unlikely]] {
                return a;
            } else {
                return m_reduce(a * pw128());
            }
        }
        modint_base(): r(0) {}
        modint_base(int64_t rr): r(rr % mod()) {
            r = std::min(r, r + mod());
            r = m_transform(r);
        }
        modint inv() const {
            return bpow(to_modint(), mod() - 2);
        }
        modint operator - () const {
            modint neg;
            neg.r = std::min(-r, mod() - r);
            return neg;
        }
        modint& operator /= (const modint &t) {
            return to_modint() *= t.inv();
        }
        modint& operator *= (const modint &t) {
            r = m_reduce(__uint128_t(r) * t.r);
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
        // Why <=> doesn't work?..
        auto operator == (const modint_base &t) const {
            return std::min(r, r - mod()) == std::min(t.r, t.r - mod());
        }
        auto operator != (const modint_base &t) const {
            return std::min(r, r - mod()) != std::min(t.r, t.r - mod());
        }
        auto operator <= (const modint_base &t) const {
            return std::min(r, r - mod()) <= std::min(t.r, t.r - mod());
        }
        auto operator >= (const modint_base &t) const {
            return std::min(r, r - mod()) >= std::min(t.r, t.r - mod());
        }
        auto operator < (const modint_base &t) const {
            return std::min(r, r - mod()) < std::min(t.r, t.r - mod());
        }
        auto operator > (const modint_base &t) const {
            return std::min(r, r - mod()) > std::min(t.r, t.r - mod());
        }
        int64_t rem() const {
            uint64_t R = getr();
            return 2 * R > (uint64_t)mod() ? R - mod() : R;
        }

        // Only use if you really know what you're doing!
        uint64_t modmod() const {return 8ULL * mod() * mod();};
        void add_unsafe(uint64_t t) {r += t;}
        void pseudonormalize() {r = std::min(r, r - modmod());}
        modint const& normalize() {
            if(r >= (uint64_t)mod()) {
                r %= mod();
            }
            return to_modint();
        }
        void setr(uint64_t rr) {r = m_transform(rr);}
        uint64_t getr() const {
            uint64_t res = m_reduce(r);
            return std::min(res, res - mod());
        }
        void setr_direct(uint64_t rr) {r = rr;}
        uint64_t getr_direct() const {return r;}
    private:
        uint64_t r;
        modint& to_modint() {return static_cast<modint&>(*this);}
        modint const& to_modint() const {return static_cast<modint const&>(*this);}
    };
    template<typename modint>
    std::istream& operator >> (std::istream &in, modint_base<modint> &x) {
        uint64_t r;
        auto &res = in >> r;
        x.setr(r);
        return res;
    }
    template<typename modint>
    std::ostream& operator << (std::ostream &out, modint_base<modint> const& x) {
        return out << x.getr();
    }

    template<typename modint>
    concept modint_type = std::is_base_of_v<modint_base<modint>, modint>;

    template<int64_t m>
    struct modint: modint_base<modint<m>> {
        static constexpr uint64_t im = m % 2 ? inv64(-m) : 0;
        static constexpr uint64_t r2 = __uint128_t(-1) % m + 1;
        static constexpr int64_t mod() {return m;}
        static constexpr uint64_t imod() {return im;}
        static constexpr __uint128_t pw128() {return r2;}
        using Base = modint_base<modint<m>>;
        using Base::Base;
    };

    struct dynamic_modint: modint_base<dynamic_modint> {
        static int64_t mod() {return m;}
        static uint64_t imod() {return im;}
        static __uint128_t pw128() {return r2;}
        static void switch_mod(int64_t nm) {
            m = nm;
            im = m % 2 ? inv64(-m) : 0;
            r2 = __uint128_t(-1) % m + 1;
        }
        using Base = modint_base<dynamic_modint>;
        using Base::Base;

        // Wrapper for temp switching
        auto static with_mod(int64_t tmp, auto callback) {
            struct scoped {
                int64_t prev = mod();
                ~scoped() {switch_mod(prev);}
            } _;
            switch_mod(tmp);
            return callback();
        }
    private:
        static int64_t m;
        static uint64_t im, r1, r2;
    };
    int64_t dynamic_modint::m = 1;
    uint64_t dynamic_modint::im = -1;
    uint64_t dynamic_modint::r2 = 0;
}
#endif // CP_ALGO_MATH_MODINT_HPP
