#ifndef CP_ALGO_ALGEBRA_NUMBER_THEORY_HPP
#define CP_ALGO_ALGEBRA_NUMBER_THEORY_HPP
#include "../random/rng.hpp"
#include "affine.hpp"
#include "modint.hpp"
#include <algorithm>
#include <optional>
#include <vector>
#include <bit>
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
    // https://en.wikipedia.org/wiki/Miller–Rabin_primality_test
    bool is_prime(uint64_t m) {
        if(m == 1 || m % 2 == 0) {
            return m == 2;
        }
        // m - 1 = 2^s * d
        int s = std::countr_zero(m - 1);
        auto d = (m - 1) >> s;
        using base = dynamic_modint;
        auto test = [&](base x) {
            x = bpow(x, d);
            if(std::abs(x.rem()) <= 1) {
                return true;
            }
            for(int i = 1; i < s && x != -1; i++) {
                x *= x;
            }
            return x == -1;
        };
        return base::with_switched_mod(m, [&](){
            // Works for all m < 2^64: https://miller-rabin.appspot.com
            return std::ranges::all_of(std::array{
                2, 325, 9375, 28178, 450775, 9780504, 1795265022
            }, test);
        });
    }
    // https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
    void factorize(uint64_t m, std::vector<int64_t> &res) {
        if(m % 2 == 0) {
            factorize(m / 2, res);
            res.push_back(2);
        } else if(is_prime(m)) {
            res.push_back(m);
        } else if(m > 1) {
            using base = dynamic_modint;
            base::with_switched_mod(m, [&]() {
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
    auto factorize(int64_t m) {
        std::vector<int64_t> res;
        factorize(m, res);
        return res;
    }
}
#endif // CP_ALGO_ALGEBRA_NUMBER_THEORY_HPP
