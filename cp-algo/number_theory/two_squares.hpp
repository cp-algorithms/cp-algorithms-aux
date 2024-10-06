#ifndef CP_ALGO_NUMBER_THEORY_TWO_SQUARES_HPP
#define CP_ALGO_NUMBER_THEORY_TWO_SQUARES_HPP
#include "euler.hpp"
#include <cassert>
#include <complex>
#include <utility>
#include <vector>
#include <map>
namespace cp_algo::math {
    using gaussint = std::complex<int64_t>;
    gaussint two_squares_prime_any(int64_t p) {
        if(p == 2) {
            return gaussint(1, 1);
        }
        assert(p % 4 == 1);
        using base = dynamic_modint;
        return base::with_mod(p, [&](){
            base g = primitive_root(p);
            int64_t i = bpow(g, (p - 1) / 4).getr();
            int64_t q0 = 1, q1 = 0;
            int64_t r = i, m = p;
            // TODO: Use library contfrac?
            do {
                int64_t d = r / m;
                q0 = std::exchange(q1, q0 + d * q1);
                r = std::exchange(m, r % m);
            } while(q1 < p / q1);
            return gaussint(q0, (base(i) * base(q0)).rem());
        });
    }

    std::vector<gaussint> two_squares_all(int64_t n) {
        if(n == 0) {
            return {0};
        }
        auto primes = factorize(n);
        std::map<int64_t, int> cnt;
        for(auto p: primes) {
            cnt[p]++;
        }
        // 1, -1, i, -i
        std::vector<gaussint> res = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
        for(auto [p, c]: cnt) {
            std::vector<gaussint> nres;
            if(p % 4 == 3) {
                if(c % 2 == 0) {
                    auto mul = bpow(gaussint(p), c / 2);
                    for(auto p: res) {
                        nres.push_back(p * mul);
                    }
                }
            } else if(p % 4 == 1) {
                gaussint base = two_squares_prime_any(p);
                for(int i = 0; i <= c; i++) {
                    auto mul = bpow(base, i) * bpow(conj(base), c - i);
                    for(auto p: res) {
                        nres.push_back(p * mul);
                    }
                }
            } else if(p % 4 == 2) {
                auto mul = bpow(gaussint(1, 1), c);
                for(auto p: res) {
                    nres.push_back(p * mul);
                }
            }
            res = nres;
        }
        std::vector<gaussint> nres;
        for(auto p: res) {
            if(p.real() >= 0 && p.imag() >= 0) {
                nres.push_back(p);
            }
        }
        return nres;
    }
}
#endif // CP_ALGO_NUMBER_THEORY_TWO_SQUARES_HPP
