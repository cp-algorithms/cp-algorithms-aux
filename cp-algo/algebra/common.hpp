#ifndef CP_ALGO_ALGEBRA_COMMON_HPP
#define CP_ALGO_ALGEBRA_COMMON_HPP
#include <functional>
#include <cstdint>
namespace cp_algo::algebra {
    const int maxn = 1 << 20;
    const int magic = 250; // threshold for sizes to run the naive algo

    auto bpow(auto x, int64_t n, auto ans, auto op) {
        for(; n; n /= 2, x = op(x, x)) {
            if(n % 2) {
                ans = op(ans, x);
            }
        }
        return ans;
    }
    auto bpow(auto x, int64_t n, auto ans) {
        return bpow(x, n, ans, std::multiplies{});
    }
    template<typename T>
    T bpow(T const& x, int64_t n) {
        return bpow(x, n, T(1));
    }

    template<typename T>
    T fact(int n) {
        static T F[maxn];
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
    
    template<typename T>
    T rfact(int n) {
        static T F[maxn];
        static bool init = false;
        if(!init) {
            F[maxn - 1] = T(1) / fact<T>(maxn - 1);
            for(int i = maxn - 2; i >= 0; i--) {
                F[i] = F[i + 1] * T(i + 1);
            }
            init = true;
        }
        return F[n];
    }

    template<typename T>
    T small_inv(int n) {
        static T F[maxn];
        static bool init = false;
        if(!init) {
            for(int i = 1; i < maxn; i++) {
                F[i] = rfact<T>(i) * fact<T>(i - 1);
            }
            init = true;
        }
        return F[n];
    }
}
#endif // CP_ALGO_ALGEBRA_COMMON_HPP
