#ifndef CP_ALGO_MATH_FACTORIZE_HPP
#define CP_ALGO_MATH_FACTORIZE_HPP
#include "primality.hpp"
#include "../random/rng.hpp"
#include <generator>
namespace cp_algo::math {
    // https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
    std::generator<uint64_t> factorize(uint64_t m) {
        if(m % 2 == 0) {
            co_yield std::ranges::elements_of(factorize(m / 2));
            co_yield 2;
        } else if(is_prime(m)) {
            co_yield m;
        } else if(m > 1) {
            using base = dynamic_modint<int64_t>;
            auto g = base::with_mod(m, [&]() {
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
                return g.getr();
            });
            co_yield std::ranges::elements_of(factorize(g));
            co_yield std::ranges::elements_of(factorize(m / g));
        }
    }
}
#endif // CP_ALGO_MATH_FACTORIZE_HPP
