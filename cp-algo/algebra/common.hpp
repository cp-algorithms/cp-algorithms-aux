#ifndef CP_ALGO_ALGEBRA_COMMON_HPP
#define CP_ALGO_ALGEBRA_COMMON_HPP
#include <functional>
#include <cstdint>
namespace cp_algo::algebra {
#ifdef CP_ALGO_MAXN
    const int maxn = CP_ALGO_MAXN;
#else
    const int maxn = 1 << 20;
#endif
    const int magic = 250; // threshold for sizes to run the naive algo

    auto bpow(auto const& x, int64_t n, auto const& one, auto op) {
        if(n == 0) {
            return one;
        } else {
            auto t = bpow(x, n / 2, one, op);
            t = op(t, t);
            if(n % 2) {
                t = op(t, x);
            }
            return t;
        }
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
    template<typename T>
    T rfact(int n) {
        static std::vector<T> F(maxn);
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
    T nCr(int n, int r) {
        if(r < 0 || r > n) {
            return T(0);
        } else {
            return fact<T>(n) * rfact<T>(r) * rfact<T>(n-r);
        }
    }
}
#endif // CP_ALGO_ALGEBRA_COMMON_HPP
