#ifndef CP_ALGO_NUMBER_THEORY_TWO_SQUARES_HPP
#define CP_ALGO_NUMBER_THEORY_TWO_SQUARES_HPP
#include "euler.hpp"
#include "../util/complex.hpp"
#include "../util/big_alloc.hpp"
#include <cassert>
#include <utility>
#include <vector>
#include <map>
namespace cp_algo::math {
    template<typename T>
    using gaussint = complex<T>;
    // Prime p == 2 or p == 1 (mod 4), with p <= max(make_signed_t<_Int>). Coefficients are signed.
    template<typename _Int>
    auto two_squares_prime_any(_Int p) {
        using Int = std::make_signed_t<_Int>;
        if(p == 2) {
            return gaussint<Int>{1, 1};
        }
        assert(p % 4 == 1);
        using base = dynamic_modint<Int>;
        return base::with_mod(p, [&](){
            base g = primitive_root(Int(p));
            int64_t i = bpow(g, (p - 1) / 4).getr();
            int64_t q0 = 1, q1 = 0;
            int64_t r = i, m = p;
            // TODO: Use library contfrac?
            do {
                int64_t d = r / m;
                q0 = std::exchange(q1, q0 + d * q1);
                r = std::exchange(m, r % m);
            } while(q1 < p / q1);
            return gaussint<Int>{Int(q0), Int((base(i) * base(q0)).rem())};
        });
    }

    // 0 <= n <= max(make_signed_t<Int>); signed Gaussian coefficients avoid wrapping conjugates.
    template<typename Int>
    big_vector<gaussint<std::make_signed_t<Int>>> two_squares_all(Int n) {
        if(n == 0) {
            return {0};
        }
        auto primes = factorize(n);
        big_map<Int, int> cnt;
        for(auto p: primes) {
            cnt[p]++;
        }
        big_vector<gaussint<std::make_signed_t<Int>>> res = {1};
        for(auto [p, c]: cnt) {
            if(p % 4 == 3 && c % 2) {
                return {};
            }
            if(p % 4 != 1) {
                auto mul = p == 2 ? bpow(gaussint<std::make_signed_t<Int>>(1, 1), c) : bpow(gaussint<std::make_signed_t<Int>>(p), c / 2);
                for(auto &value: res) {
                    value *= mul;
                }
            } else {
                big_vector<gaussint<std::make_signed_t<Int>>> nres;
                auto base = two_squares_prime_any(p);
                for(int i = 0; i <= c; i++) {
                    auto mul = bpow(base, i) * bpow(conj(base), c - i);
                    for(auto value: res) {
                        nres.push_back(value * mul);
                    }
                }
                res = std::move(nres);
            }
        }
        big_vector<gaussint<std::make_signed_t<Int>>> nres;
        for(auto p: res) {
            while(p.real() < 0 || p.imag() < 0) {
                p *= gaussint<std::make_signed_t<Int>>(0, 1);
            }
            nres.push_back(p);
            if(!p.real() || !p.imag()) {
                nres.emplace_back(p.imag(), p.real());
            }
        }
        return nres;
    }
}
#endif // CP_ALGO_NUMBER_THEORY_TWO_SQUARES_HPP
