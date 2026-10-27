// @brief Dense Convolution
#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"
#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/io.hpp"
#include "cp-algo/math/fft.hpp"

using cp_algo::io::read_reduced;
using cp_algo::io::write;

using namespace cp_algo::math;
using base = modint<998244353>;

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    int n, m;
    std::cin >> n >> m;
    cp_algo::big_vector<base> a(n), b(m);
    read_reduced(a);
    read_reduced(b);
    fft::mul(a, b);
    write(a);
    std::cout << '\n';
}
