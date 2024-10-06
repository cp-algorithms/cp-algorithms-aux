#ifndef CP_ALGO_MATH_NUMBER_THEORY_HPP
#define CP_ALGO_MATH_NUMBER_THEORY_HPP
#include "../random/rng.hpp"
#include "affine.hpp"
#include "factorize.hpp"
#include <algorithm>
#include <optional>
#include <vector>
#include <bit>
namespace cp_algo::math {
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
