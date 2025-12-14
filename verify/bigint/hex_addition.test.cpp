// @brief Addition of Hex Big Integers
#define PROBLEM "https://judge.yosupo.jp/problem/addition_of_hex_big_integers"
#pragma GCC optimize("O3,unroll-loops")
#include <bits/allocator.h>
#pragma GCC target("avx2")
#define CP_ALGO_CHECKPOINT
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/bigint.hpp"
#include "cp-algo/util/checkpoint.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::math;

void solve() {
    bigint<x16> a, b;
    cin >> a >> b;
    cout << (a += b) << "\n";
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    cin >> t;
    while(t--) {
        solve();
    }
    cp_algo::checkpoint<1>();
}
