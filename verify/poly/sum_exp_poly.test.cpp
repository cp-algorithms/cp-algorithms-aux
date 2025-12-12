// @brief $\sum\limits_{i=0}^\infty r^i i^d$
#define PROBLEM "https://judge.yosupo.jp/problem/sum_of_exponential_times_polynomial_limit"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#define CP_ALGO_MAXN 1 << 24
#include "cp-algo/math/poly.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
typedef modint<mod> base;
typedef poly_t<base> polyn;

polyn xm1k(size_t k) {
    polyn::Vector ans(k+1);
    for(size_t i = 0; i <= k; i++) {
        ans[i] = (k - i) & 1 ? -binom<base>(k, i) : binom<base>(k, i);
    }
    return ans;
}

void solve() {
    int r, d;
    cin >> r >> d;
    polyn num = polyn(bpow(base(1-r), d+1)) - bpow(base(r), d+1) * xm1k(d+1);
    polyn den = polyn(bpow(base(1-r), d+2)) - bpow(base(1-r), d+1) * r * (polyn::xk(1) - polyn(base(1)));
    polyn H = num / den;
    base ans = 0;

    vector<base> id(d+1);
    vector<int> lp(d+1);
    id[0] = bpow(base(0), d);
    id[1] = 1;
    for(int i = 2; i <= d; i++) {
        if(lp[i] == 0) {
            id[i] = bpow(base(i), d);
            for(int j = i; j <= d; j += i) {
                lp[j] = i;
            }
        } else {
            id[i] = id[lp[i]] * id[i / lp[i]];
        }
    }

    for(int i = 0; i <= d; i++) {
        ans += H[i] * id[i];
    }
    cout << ans << endl;
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t;
    t = 1;// cin >> t;
    while(t--) {
        solve();
    }
}
