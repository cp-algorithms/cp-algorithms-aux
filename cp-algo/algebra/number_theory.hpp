#ifndef CP_ALGO_ALGEBRA_NUMBER_THEORY_HPP
#define CP_ALGO_ALGEBRA_NUMBER_THEORY_HPP
#include "../random/rng.hpp"
#include "affine.hpp"
#include "modint.hpp"
#include <optional>
namespace cp_algo::algebra {
    // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm
    template<typename base>
    requires(std::is_base_of_v<modint_base<base>, base>)
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
#endif // CP_ALGO_ALGEBRA_NUMBER_THEORY_HPP
