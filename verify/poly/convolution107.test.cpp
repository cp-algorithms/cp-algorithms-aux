// @brief Convolution mod $10^9+7$
#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod_1000000007"
#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#define CP_ALGO_CHECKPOINT
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/io.hpp"
#include "cp-algo/math/fft.hpp"

using cp_algo::io::read_reduced;

using namespace std;
using namespace cp_algo::math;

const int mod = 1e9 + 7;
using base = modint<mod>;

void solve() {
    int n, m;
    cin >> n >> m;
    cp_algo::big_vector<base> a(n), b(m);
    read_reduced(a);
    read_reduced(b);
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
