// @brief Lazy FPS: Product with a Known Polynomial
#define PROBLEM "https://judge.yosupo.jp/problem/convolution_mod"
#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/io.hpp"
#include "cp-algo/math/fps.hpp"
using cp_algo::io::read_reduced;
using namespace cp_algo::math;
using base = modint<998244353>;
using polyn = poly_t<base>;
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    int n, m;
    std::cin >> n >> m;
    polyn::Vector a(n), b(m);
    read_reduced(a);
    read_reduced(b);
    fps<base> stream([b = std::move(b)](size_t i, auto const&) {
        return i < b.size() ? b[i] : base(0);
    });
    (polyn(std::move(a)) * stream).prefix(n + m - 1).print(n + m - 1);
}
