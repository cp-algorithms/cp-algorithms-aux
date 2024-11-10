#ifndef CP_ALGO_MATH_FACTORIZE_HPP
#define CP_ALGO_MATH_FACTORIZE_HPP
#include "primality.hpp"
#include "../random/rng.hpp"
namespace cp_algo::math {
    // https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
    std::basic_string<uint64_t> factorize(uint64_t m) {
        if(m % 2 == 0) {
            return factorize(m / 2) + (uint64_t)2;
        } else if(is_prime(m)) {
            return {m};
        } else if(m > 1) {
            using base = dynamic_modint<int64_t>;
            return base::with_mod(m, [&]() {
                base t = random::rng();
                auto f = [&](auto x) {
                    return x * x + t;
                };
                base x, y;
                base g = 1;
                while(g == 1) {
                    for(int i = 0; i < 64; i++) {
                        x = f(x);
                        y = f(f(y));
                        if(x == y) [[unlikely]] {
                            t = random::rng();
                            x = y = 0;
                        } else {
                            base t = g * (x - y);
                            g = t == 0 ? g : t;
                        }
                    }
                    g = std::gcd(g.getr(), m);
                }
                return factorize(g.getr()) + factorize(m / g.getr());
            });
        } else {
            return {};
        }
    }
}
#endif // CP_ALGO_MATH_FACTORIZE_HPP
