// @brief Dense Convolution
#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"
#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/fft.hpp"

using namespace cp_algo::math;
using base = modint<998244353>;

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    int n, m;
    std::cin >> n >> m;
    cp_algo::big_vector<base> a(n), b(m);
    for(auto &x: a) {std::cin >> x;}
    for(auto &x: b) {std::cin >> x;}
    fft::mul(a, b);
    for(auto x: a) {std::cout << x << ' ';}
    std::cout << '\n';
}
