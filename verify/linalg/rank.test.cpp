// @brief Rank of Matrix
#define PROBLEM "https://judge.yosupo.jp/problem/matrix_rank"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <bits/stdc++.h>
//#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/io.hpp"
#include "cp-algo/linalg/matrix.hpp"

using cp_algo::io::read_reduced;

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
        read_reduced(A);
    } else {
        A = matrix<base>(m, n);
        for(size_t j = 0; j < A.m(); j++) {
            for(size_t i = 0; i < A.n(); i++) {
                read_reduced(A[i][j]);
            }
        }
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
