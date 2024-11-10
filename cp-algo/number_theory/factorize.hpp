#ifndef CP_ALGO_MATH_FACTORIZE_HPP
#define CP_ALGO_MATH_FACTORIZE_HPP
#include "primality.hpp"
#include "../random/rng.hpp"
namespace cp_algo::math {
    // https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
    void factorize(uint64_t m, std::vector<int64_t> &res) {
        if(m % 2 == 0) {
            factorize(m / 2, res);
            res.push_back(2);
        } else if(is_prime(m)) {
            res.push_back(m);
        } else if(m > 1) {
            using base = dynamic_modint<>;
            base::with_mod(m, [&]() {
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
                factorize(g.getr(), res);
                factorize(m / g.getr(), res);
            });
        }
    }
    std::vector<int64_t> factorize(int64_t m) {
        std::vector<int64_t> res;
        factorize(m, res);
        return res;
    }
}
#endif // CP_ALGO_MATH_FACTORIZE_HPP
