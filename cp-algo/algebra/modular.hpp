#ifndef CP_ALGO_ALGEBRA_MODULAR_HPP
#define CP_ALGO_ALGEBRA_MODULAR_HPP
#include "../random/rng.hpp"
#include "affine.hpp"
#include "common.hpp"
#include <algorithm>
#include <iostream>
#include <optional>
namespace cp_algo::algebra {
    template<int m>
    struct modular {
        // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm
        std::optional<modular> sqrt() const {
            if(r == 0) {
                return 0;
            } else if(bpow(*this, (m - 1) / 2) != modular(1)) {
                return std::nullopt;
            } else {
                while(true) {
                    modular z = random::rng();
                    if(z * z == *this) {
                        return z;
                    }
                    lin<modular> x(1, z, *this); // x + z (mod x^2 - b)
                    x = bpow(x, (m - 1) / 2, lin<modular>(0, 1, *this));
                    if(x.a != modular(0)) {
                        return x.a.inv();
                    }
                }
            }
        }
        
        uint64_t r;
        constexpr modular(): r(0) {}
        constexpr modular(int64_t rr): r(rr % m) {r = std::min<uint64_t>(r, r + m);}
        modular inv() const {return bpow(*this, m - 2);}
        modular operator - () const {return std::min(-r, m - r);}
        modular operator * (const modular &t) const {return r * t.r;}
        modular operator / (const modular &t) const {return *this * t.inv();}
        modular& operator += (const modular &t) {r += t.r; r = std::min<uint64_t>(r, r - m); return *this;}
        modular& operator -= (const modular &t) {r -= t.r; r = std::min<uint64_t>(r, r + m); return *this;}
        modular operator + (const modular &t) const {return modular(*this) += t;}
        modular operator - (const modular &t) const {return modular(*this) -= t;}
        modular& operator *= (const modular &t) {return *this = *this * t;}
        modular& operator /= (const modular &t) {return *this = *this / t;}
        
        auto operator <=> (const modular &t) const = default;
        
        explicit operator int() const {return r;}
        int64_t rem() const {return 2 * r > m ? r - m : r;}

        static constexpr uint64_t mm = (uint64_t)m * m;
        void add_unsafe(uint64_t t) {r += t; r = std::min<uint64_t>(r, r - mm);}
        modular& normalize() {if(r >= m) r %= m; return *this;}
    };
    
    template<int m>
    std::istream& operator >> (std::istream &in, modular<m> &x) {
        return in >> x.r;
    }
    
    template<int m>
    std::ostream& operator << (std::ostream &out, modular<m> const& x) {
        return out << x.r % m;
    }
}
#endif // CP_ALGO_ALGEBRA_MODULAR_HPP
