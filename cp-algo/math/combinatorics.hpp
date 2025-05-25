#ifndef CP_ALGO_MATH_COMBINATORICS_HPP
#define CP_ALGO_MATH_COMBINATORICS_HPP
#include "common.hpp"
#include <cassert>
namespace cp_algo::math {
    // fact/rfact/small_inv are caching
    // Beware of usage with dynamic mod
    template<typename T>
    T fact(auto n) {
        static std::vector<T> F(maxn);
        static bool init = false;
        if(!init) {
            F[0] = T(1);
            for(int i = 1; i < maxn; i++) {
                F[i] = F[i - 1] * T(i);
            }
            init = true;
        }
        return F[n];
    }
    // Only works for modint types
    template<typename T>
    T rfact(auto n) {
        static std::vector<T> F(maxn);
        static bool init = false;
        if(!init) {
            int t = (int)std::min<int64_t>(T::mod(), maxn) - 1;
            F[t] = T(1) / fact<T>(t);
            for(int i = t - 1; i >= 0; i--) {
                F[i] = F[i + 1] * T(i + 1);
            }
            init = true;
        }
        return F[n];
    }
    template<typename T>
    T small_inv(auto n) {
        static std::vector<T> F(maxn);
        static bool init = false;
        if(!init) {
            for(int i = 1; i < maxn; i++) {
                F[i] = rfact<T>(i) * fact<T>(i - 1);
            }
            init = true;
        }
        return F[n];
    }
    template<typename T>
    T binom_large(T n, auto r) {
        assert(r < maxn);
        T ans = 1;
        for(decltype(r) i = 0; i < r; i++) {
            ans = ans * T(n - i) * small_inv<T>(i + 1);
        }
        return ans;
    }
    template<typename T>
    T binom(auto n, auto r) {
        if(r < 0 || r > n) {
            return T(0);
        } else if(n >= maxn) {
            return binom_large(T(n), r);
        } else {
            return fact<T>(n) * rfact<T>(r) * rfact<T>(n - r);
        }
    }
}
#endif // CP_ALGO_MATH_COMBINATORICS_HPP
