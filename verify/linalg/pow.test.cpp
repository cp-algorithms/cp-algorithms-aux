// @brief Pow of Matrix
#define PROBLEM "https://judge.yosupo.jp/problem/pow_of_matrix"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("tune=native")
#include "cp-algo/linalg/matrix.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math::linalg;
using namespace cp_algo::math;

const int mod = 998244353;

void solve() {
    int n;
    uint64_t k;
    cin >> n >> k;
    matrix<modint<mod>> a(n, n);
    a.read();
    a.pow(k).print();
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