#ifndef CP_ALGO_MATH_NUMBER_THEORY_HPP
#define CP_ALGO_MATH_NUMBER_THEORY_HPP
#include "../random/rng.hpp"
#include "affine.hpp"
#include "modint.hpp"
#include <algorithm>
#include <optional>
#include <vector>
#include <bit>
namespace cp_algo::math {
    std::vector<int64_t> factorize(int64_t m);

    int64_t euler_phi(int64_t m) {
        auto primes = factorize(m);
        std::ranges::sort(primes);
        auto [from, to] = std::ranges::unique(primes);
        primes.erase(from, to);
        int64_t ans = m;
        for(auto it: primes) {
            ans -= ans / it;
        }
        return ans;
    }
    template<modint_type base>
    int64_t period(base x) {
        auto ans = euler_phi(base::mod());
        base x0 = bpow(x, ans);
        for(auto t: factorize(ans)) {
            while(ans % t == 0 && x0 * bpow(x, ans / t) == x0) {
                ans /= t;
            }
        }
        return ans;
    }
    // Find min non-negative x s.t. a*b^x = c (mod m)
    std::optional<uint64_t> discrete_log(int64_t b, int64_t c, uint64_t m, int64_t a = 1) {
        if(std::abs(a - c) % m == 0) {
            return 0;
        }
        if(std::gcd(a, m) != std::gcd(a * b, m)) {
            auto res = discrete_log(b, c, m, a * b % m);
            return res ? std::optional(*res + 1) : res;
        }
        // a * b^x is periodic here
        using base = dynamic_modint;
        return base::with_mod(m, [&]() -> std::optional<uint64_t> {
            size_t sqrtmod = std::max<size_t>(1, std::sqrt(m) / 2);
            std::unordered_map<int64_t, int> small;
            base cur = a;
            for(size_t i = 0; i < sqrtmod; i++) {
                small[cur.getr()] = i;
                cur *= b;
            }
            base step = bpow(base(b), sqrtmod);
            cur = 1;
            for(size_t k = 0; k < m; k += sqrtmod) {
                auto it = small.find((base(c) * cur).getr());
                if(it != end(small)) {
                    auto cand = base::with_mod(period(base(b)), [&](){
                        return base(it->second - k);
                    }).getr();
                    if(base(a) * bpow(base(b), cand) == base(c)) {
                        return cand;
                    } else {
                        return std::nullopt;
                    }
                }
                cur *= step;
            }
            return std::nullopt;
        });
    }
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
        return base::with_mod(m, [&](){
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
    int64_t primitive_root(int64_t p) {
        using base = dynamic_modint;
        return base::with_mod(p, [p](){
            base t = 1;
            while(period(t) != p - 1) {
                t = random::rng();
            }
            return t.getr();
        });
    }
}
#endif // CP_ALGO_MATH_NUMBER_THEORY_HPP
