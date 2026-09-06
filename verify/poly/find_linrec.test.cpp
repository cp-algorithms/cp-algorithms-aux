// @brief Find Linear Recurrence
#define PROBLEM "https://judge.yosupo.jp/problem/find_linear_recurrence"
#pragma GCC optimize("O3,unroll-loops")
//#include <bits/allocator.h>
#pragma GCC target("avx2")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/poly/recurrence.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n;
    cin >> n;
    polyn::Vector a(n);
    for(auto &it: a) {cin >> it;}
    auto Q = min_rec(polyn(std::move(a)), n);
    int d = Q.deg();
    cout << d << endl;
    (-Q / Q[d]).reverse().div_xk(1).print(d);
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
