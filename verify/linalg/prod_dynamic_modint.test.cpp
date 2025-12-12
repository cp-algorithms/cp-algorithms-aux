// @brief Matrix Product (dynamic modint)
#define PROBLEM "https://judge.yosupo.jp/problem/matrix_product"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/linalg/matrix.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::linalg;
using namespace cp_algo::math;
using base = dynamic_modint<int64_t>;

const int64_t mod = 998244353;

void solve() {
    base::switch_mod(mod);
    int n, m, k;
    cin >> n >> m >> k;
    matrix<base> a(n, m), b(m, k);
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