namespace algebra { // common
    const int maxn = 1 << 20;
    const int magic = 250; // threshold for sizes to run the naive algo
    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count()); 

    auto bpow(auto x, int64_t n, auto ans) {
        for(; n; n /= 2, x = x * x) {
            if(n % 2) {
                ans = ans * x;
            }
        }
        return ans;
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
