// @brief System of Linear Equations
#define PROBLEM "https://judge.yosupo.jp/problem/system_of_linear_equations"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("avx2,tune=native")
#include "cp-algo/algebra/matrix.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::algebra;

const int mod = 998244353;

void solve() {
    int n, m;
    cin >> n >> m;
    matrix<modular<mod>> A(n, m), b(n, 1);
    A.read();
    b.read();
    auto x = A.solve(b);
    if(!x) {
        cout << -1 << "\n";
    } else {
        auto [sol, basis] = *x;
        cout << basis.n << "\n";
        sol.print();
        basis.print();
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