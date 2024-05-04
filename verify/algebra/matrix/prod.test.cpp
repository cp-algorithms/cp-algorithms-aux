// @brief Matrix Product
#define PROBLEM "https://judge.yosupo.jp/problem/matrix_product"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("avx2,tune=native")
#include "cp-algo/linalg/matrix.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::linalg;
using namespace cp_algo::algebra;

const int mod = 998244353;

void solve() {
    int n, m, k;
    cin >> n >> m >> k;
    matrix<modint<mod>> a(n, m), b(m, k);
    a.read();
    b.read();
    (a * b).print();
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