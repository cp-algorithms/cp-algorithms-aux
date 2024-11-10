// @brief Adjugate Matrix
#define PROBLEM "https://judge.yosupo.jp/problem/adjugate_matrix"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/linalg/matrix.hpp"
#include <bits/stdc++.h>

const int64_t mod = 998244353;

using namespace std;
using cp_algo::math::modint;
using cp_algo::linalg::matrix;


void solve() {
    int n;
    cin >> n;
    matrix<modint<mod>> A(n + 1);
    for(int i: views::iota(0, n)) {
        for(int j: views::iota(0, n)) {
            cin >> A[i][j];
        }
    }
    for(int i: views::iota(0, n)) {
        A[i][n] = cp_algo::random::rng();
        A[n][i] = cp_algo::random::rng();
    }
    auto [D, Ai] = A.inv();
    for(int i: views::iota(0, n)) {
        for(int j: views::iota(0, n)) {
            if(D != 0) {
                auto res = Ai[n][n] * Ai[i][j] - Ai[i][n] * Ai[n][j];
                cout << res * D << " \n"[j + 1 == n];
            } else {
                cout << 0 << " \n"[j + 1 == n];
            }
            
        }
    }
}

int main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    solve();
    return 0;
}
