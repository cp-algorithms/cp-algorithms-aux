// @brief Pow of Power Series (Sparse)
#define PROBLEM "https://judge.yosupo.jp/problem/pow_of_formal_power_series_sparse"
#define CP_ALGO_MAXN (1 << 20)
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/poly/sparse/pow.hpp"
using namespace std;
using namespace cp_algo::math;
using base = modint<998244353>;
using polyn = poly_t<base>;
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n, k;
    cin >> n >> k;
    int64_t m;
    cin >> m;
    polyn::Vector a(n);
    for(int j = 0; j < k; j++) {
        int i;
        cin >> i;
        cin >> a[i];
    }
    pow_sparse(polyn(std::move(a)), m, n).print(n);
}
