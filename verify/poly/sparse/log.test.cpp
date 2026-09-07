// @brief Log of Power Series (Sparse)
#define PROBLEM "https://judge.yosupo.jp/problem/log_of_formal_power_series_sparse"
#define CP_ALGO_MAXN (1 << 20)
#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/poly/sparse/log.hpp"
using namespace std;
using namespace cp_algo::math;
using base = modint<998244353>;
using polyn = poly_t<base>;
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n, k;
    cin >> n >> k;
    polyn::Vector a(n);
    for(int j = 0; j < k; j++) {
        int i;
        cin >> i;
        cin >> a[i];
    }
    log_sparse(polyn(std::move(a)), n).print(n);
}
