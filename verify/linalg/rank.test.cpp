// @brief Rank of Matrix
#define PROBLEM "https://judge.yosupo.jp/problem/matrix_rank"
#pragma GCC optimize("Ofast,unroll-loops")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/linalg/matrix.hpp"

using namespace std;
using namespace cp_algo::math;
using namespace cp_algo::linalg;

const int64_t mod = 998244353;
using base = modint<mod>;

void solve() {
    int n, m;
    cin >> n >> m;
    matrix<base> A;
    if(n < m) {
        A = matrix<base>(n, m);
        A.read();
    } else {
        A = matrix<base>(m, n);
        A.read_transposed();
    }
    cout << A.rank() << "\n";
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
