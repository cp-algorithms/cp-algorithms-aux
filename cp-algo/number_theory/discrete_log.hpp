#ifndef CP_ALGO_NUMBER_THEORY_DISCRETE_LOG_HPP
#define CP_ALGO_NUMBER_THEORY_DISCRETE_LOG_HPP
#include "euler.hpp"
namespace cp_algo::math {
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
                    auto cand = base::with_mod(period(base(b)), [&]() {
                        return base(it->second - k).getr();
                    });
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
}
#endif // CP_ALGO_NUMBER_THEORY_DISCRETE_LOG_HPP
