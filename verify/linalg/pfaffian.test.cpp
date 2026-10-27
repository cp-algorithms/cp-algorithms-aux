// @brief Pfaffian of Matrix
#define PROBLEM "https://judge.yosupo.jp/problem/pfaffian_of_matrix"
#include <bits/stdc++.h>
#include "cp-algo/util/io.hpp"
#include "cp-algo/linalg/matrix.hpp"

using cp_algo::io::read_reduced;

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    size_t n;
    std::cin >> n;
    using base = cp_algo::math::modint<998244353LL>;
    cp_algo::linalg::matrix<base> a(2 * n);
    read_reduced(a);
    std::cout << a.pfaffian() << '\n';
}
