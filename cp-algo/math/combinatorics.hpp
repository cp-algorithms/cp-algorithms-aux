#ifndef CP_ALGO_MATH_COMBINATORICS_HPP
#define CP_ALGO_MATH_COMBINATORICS_HPP
#include "common.hpp"
namespace cp_algo::math {
    template<typename T>
    T binom_large(T n, int r) {
        assert(r < math::maxn);
        T ans = 1;
        for(int i = 0; i < r; i++) {
            ans = ans * T(n - i) * small_inv<T>(i + 1);
        }
        return ans;
    }
    template<typename T>
    T binom(int n, int r) {
        if(r < 0 || r > n) {
            return T(0);
        } else if(n > math::maxn) {
            return binom_large(n, r);
        } else {
            return math::fact<T>(n) * math::rfact<T>(r) * math::rfact<T>(n-r);
        }
    }
}
#endif // CP_ALGO_MATH_COMBINATORICS_HPP
