// @brief Enumerate Primes
#define PROBLEM "https://judge.yosupo.jp/problem/enumerate_primes"
#pragma GCC optimize("O3,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/checkpoint.hpp"
#include "cp-algo/math/sieve.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

void solve() {
    uint32_t N, A, B;
    cin >> N >> A >> B;
    auto primes = odd_sieve(N);
    auto cnt = count(primes) + (N >= 2);
    cp_algo::checkpoint("count");
    auto X = cnt < B ? 0 : (cnt - B + A - 1) / A;
    cout << cnt << ' ' << X << endl;
    if (X) {
        if (B == 0) {
            cout << 2 << ' ';
        }
        B = (B - 1 + A) % A;
    }
    for(size_t i = skip(primes, 0, B); i < primes.size(); i = skip(primes, i, A)) {
        cout << to_val(i) << ' ';
    }
    cout << "\n";
    cp_algo::checkpoint("print");
    cp_algo::checkpoint<1>();
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    //cin >> t;
    while(t--) {
        solve();
    }
}
