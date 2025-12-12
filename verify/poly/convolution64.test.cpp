// @brief Convolution (Mod $2^{64}$)
#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod_2_64"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#define CP_ALGO_CHECKPOINT
#include <bits/stdc++.h>
#include "cp-algo/math/fft64.hpp"
//#include "blazingio/blazingio.min.hpp"

using namespace std;

void solve() {
    int n, m;
    cin >> n >> m;
    cp_algo::big_vector<uint64_t> a(n), b(m);
    for(auto &x : a) cin >> x;
    for(auto &x : b) cin >> x;
    cp_algo::checkpoint("read");
    cp_algo::math::fft::conv64(a, b);
    for(auto x: a) {
        cout << uint64_t(x) << " ";
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