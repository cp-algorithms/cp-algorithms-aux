#ifndef CP_ALGO_ALGEBRA_NUMBER_THEORY_HPP
#define CP_ALGO_ALGEBRA_NUMBER_THEORY_HPP
#include "../random/rng.hpp"
#include "affine.hpp"
#include "modint.hpp"
#include <optional>
namespace cp_algo::algebra {
    // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm
    template<modint_type base>
    std::optional<base> sqrt(base b) {
        if(b == base(0)) {
            return base(0);
        } else if(bpow(b, (b.mod() - 1) / 2) != base(1)) {
            return std::nullopt;
        } else {
            while(true) {
                base z = random::rng();
                if(z * z == b) {
                    return z;
                }
                lin<base> x(1, z, b); // x + z (mod x^2 - b)
                x = bpow(x, (b.mod() - 1) / 2, lin<base>(0, 1, b));
                if(x.a != base(0)) {
                    return x.a.inv();
                }
            }
        }
    }

    template<modint_type base>
    bool is_prime_mod() {
        auto m = base::mod();
        if(m == 1 || m % 2 == 0) {
            return m == 2;
        }
        auto m1 = m - 1;
        int d = 0;
        while(m1 % 2 == 0) {
            m1 /= 2;
            d++;
        }
        auto test = [&](auto x) {
            x = bpow(x, m1);
            if(x == 0 || x == 1 || x == -1) {
                return true;
            }
            for(int i = 0; i <= d; i++) {
                if(x == -1) {
                    return true;
                }
                x *= x;
            }
            return false;
        };
        for(base b: {2, 325, 9375, 28178, 450775, 9780504, 1795265022}) {
            if(!test(b)) {
                return false;
            }
        }
        return true;
    }
    bool is_prime(int64_t m) {
        return dynamic_modint::with_switched_mod(m, [](){
            return is_prime_mod<dynamic_modint>();
        });
    }
}
#endif // CP_ALGO_ALGEBRA_NUMBER_THEORY_HPP
