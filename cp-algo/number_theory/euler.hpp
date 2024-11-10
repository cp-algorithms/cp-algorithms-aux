#ifndef CP_ALGO_NUMBER_THEORY_EULER_HPP
#define CP_ALGO_NUMBER_THEORY_EULER_HPP
#include "factorize.hpp"
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
    int64_t primitive_root(int64_t p) {
        using base = dynamic_modint<>;
        return base::with_mod(p, [p](){
            base t = 1;
            while(period(t) != p - 1) {
                t = random::rng();
            }
            return t.getr();
        });
    }
}
#endif // CP_ALGO_NUMBER_THEORY_EULER_HPP
