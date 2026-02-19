// @brief Inv of Bivariate Formal Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/inv_of_formal_power_series_2d"
#pragma GCC optimize("O3,unroll-loops")
#include <bits/allocator.h>
#define CP_ALGO_CHECKPOINT
#pragma GCC target("avx2")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/multivar_inv.hpp"

using namespace std;
using namespace cp_algo::math::fft;

const int mod = 998244353;
using base = cp_algo::math::modint<mod>;

void solve() {
    size_t n, m;
    cin >> n >> m;
    std::array<size_t, 2> dim{m, n};
    multivar<base> a(dim);
    a.read();
    auto b = multivar_inv(a, dim);
    b.print();
    cp_algo::checkpoint<1>();
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    solve();
}
