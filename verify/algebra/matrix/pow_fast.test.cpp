// @brief Pow of Matrix (Frobenius)
#define PROBLEM "https://judge.yosupo.jp/problem/pow_of_matrix"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("avx2,tune=native")
#include "cp-algo/linalg/frobenius.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::algebra;
using namespace cp_algo::linalg;

const int mod = 998244353;
using base = modular<mod>;
using polyn = poly_t<base>;

template<typename base>
auto frobenius_pow(matrix<base> A, uint64_t k) {
    using polyn = poly_t<base>;
    auto [T, Tinv, charps] = frobenius_form<full>(A);
    vector<matrix<base>> blocks;
    for(auto charp: charps) {
        matrix<base> block(charp.deg());
        auto xk = polyn::xk(1).powmod(k, charp);
        for(size_t i = 0; i < block.n(); i++) {
            ranges::copy(xk.a, begin(block[i]));
            xk = xk.mul_xk(1) % charp;
        }
        blocks.push_back(block);
    }
    auto S = matrix<base>::block_diagonal(blocks);
    return Tinv * S * T;
}

void solve() {
    size_t n;
    uint64_t k;
    cin >> n >> k;
    matrix<base> A(n);
    A.read();
    frobenius_pow(A, k).print();
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    //cin >> t;
    while(t--) {
        solve();
    }
}
