// @brief Convolution mod $10^9+7$
#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod_1000000007"
#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/fft.hpp"

using namespace std;
using namespace cp_algo::math;

const int mod = 1e9 + 7;
using base = modint<mod>;

void solve() {
    int n, m;
    cin >> n >> m;
    vector<base, cp_algo::big_alloc<base>> a(n), b(m);
    for(auto &x: a) {cin >> x;}
    for(auto &x: b) {cin >> x;}
    cp_algo::checkpoint("read");
    fft::mul(a, b);
    for(auto x: a) {cout << x << " ";}
    cp_algo::checkpoint("write");
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    while(t--) {
        solve();
    }
}
