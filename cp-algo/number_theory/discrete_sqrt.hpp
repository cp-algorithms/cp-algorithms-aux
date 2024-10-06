#ifndef CP_ALGO_NUMBER_THEORY_DISCRETE_SQRT_HPP
#define CP_ALGO_NUMBER_THEORY_DISCRETE_SQRT_HPP
#include "modint.hpp"
#include "cp-algo/random/rng.hpp"
#include "cp-algo/math/affine.hpp"
namespace cp_algo::math {
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
}
#endif // CP_ALGO_NUMBER_THEORY_SQRT_HPP
