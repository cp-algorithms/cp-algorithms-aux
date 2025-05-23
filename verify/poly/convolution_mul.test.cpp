// @brief Convolution on the Multiplicative Monoid of $\mathbb Z/p\mathbb{Z}$
#define PROBLEM "https://judge.yosupo.jp/problem/mul_modp_convolution"
#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/number_theory/euler.hpp"
#include "cp-algo/math/fft.hpp"

using namespace std;

using base = cp_algo::math::modint<998244353>;

void solve() {
    int p;
    cin >> p;
    auto g = cp_algo::math::primitive_root(p);
    vector<int> lg(p);
    int64_t cur = 1;
    for(int i = 0; i < p - 1; i++) {
        lg[cur] = i;
        cur *= g;
        cur %= p;
    }
    cp_algo::checkpoint("find lg");
    base a0, b0, as = 0, bs = 0;
    vector<base> a(p-1), b(p-1);
    cin >> a0;
    for(int i = 1; i <= p - 1; i++) {
        cin >> a[lg[i]];
        as += a[lg[i]];
    }
    cin >> b0;
    for(int i = 1; i <= p - 1; i++) {
        cin >> b[lg[i]];
        bs += b[lg[i]];
    }
    cp_algo::checkpoint("read");
    base c0 = (a0 + as) * (b0 + bs) - as * bs;
    cout << c0 << " ";
    cp_algo::math::fft::mul(a, b);
    for(size_t i = p-1; i < size(a); i++) {
        a[i - (p-1)] += a[i];
    }
    for(int i = 1; i <= p - 1; i++) {
        cout << a[lg[i]] << " ";
    }
    cp_algo::checkpoint("write");
    cp_algo::checkpoint<1>();
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    solve();
}
