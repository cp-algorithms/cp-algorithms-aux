// @brief Inv of Polynomials
#define PROBLEM "https://judge.yosupo.jp/problem/inv_of_polynomials"
#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
//#include <bits/allocator.h>
#pragma GCC target("avx2")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/poly/euclid.hpp"

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n, m;
    cin >> n >> m;
    polyn::Vector a(n), b(m);
    for(auto &it: a) {cin >> it;}
    for(auto &it: b) {cin >> it;}
    auto res = inv_mod(polyn(std::move(a)), polyn(b));
    if(res) {
        cout << res->deg() + 1 << endl;
        res->print();
    } else {
        cout << -1 << endl;
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
