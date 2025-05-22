// @brief Multidimensional Convolution (Truncated)
#define PROBLEM "https://judge.yosupo.jp/problem/multivariate_convolution"
#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/multivar.hpp"

using namespace std;
using namespace cp_algo::math::fft;

const int mod = 998244353;
using base = cp_algo::math::modint<mod>;

void solve() {
    int k;
    cin >> k;
    vector<size_t> N(k);
    for(auto &n: N) cin >> n;
    multivar<base> a(N), b(N);
    a.read();
    b.read();
    a.mul(b);
    a.print();
    cp_algo::checkpoint<1>();
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    solve();
}
