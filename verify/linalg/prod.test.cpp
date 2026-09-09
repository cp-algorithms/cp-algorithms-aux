// @brief Matrix Product
#define PROBLEM "https://judge.yosupo.jp/problem/matrix_product"
#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
//#include "blazingio/blazingio.min.hpp"
#include "cp-algo/linalg/matrix.hpp"

using namespace std;
using namespace cp_algo::linalg;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;

void solve() {
    int n, m, k;
    cin >> n >> m >> k;
    matrix<base> a(n, m), b(m, k);
    a.read();
    b.read();
    auto c = a * b;
    for(size_t i = 0; i < c.n(); i++) {
        c[i].print();
        if(i % 64 == 63) cout.flush(); // Reuse the output buffer for large matrices.
    }
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
