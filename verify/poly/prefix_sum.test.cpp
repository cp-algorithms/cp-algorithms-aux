// @brief Prefix Sum of Polynomial
#define PROBLEM "https://judge.yosupo.jp/problem/prefix_sum_of_polynomial"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#define CP_ALGO_MAXN ((1 << 19) + 1)
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/poly/transform.hpp"

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n;
    cin >> n;
    polyn::Vector f(n);
    for(auto &x : f) cin >> x;
    prefix_sum(polyn(std::move(f))).print(n + 1);

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
