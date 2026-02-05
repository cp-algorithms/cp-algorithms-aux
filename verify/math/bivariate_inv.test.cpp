// @brief Inv of Bivariate Formal Power Series
#define PROBLEM "https://judge.yosupo.jp/problem/inv_of_formal_power_series_2d"
#pragma GCC optimize("O3,unroll-loops")
#include <bits/allocator.h>
#define CP_ALGO_CHECKPOINT
#pragma GCC target("avx2")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/multivar.hpp"

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
    size_t degree = 1;
    multivar<base> b(std::array{1, 1});
    b.data[0] = a.data[0].inv();
    while(degree < n + m) {
        degree *= 2;
        std::array next_dim{std::min(m, degree), std::min(n, degree)};

        multivar<base> c = a.truncated(next_dim);
        b.extend(next_dim);
        c.mul(b);
        for(auto &x: c.data) {
            x = -x;
        }
        c.data[0] += base(2);
        b.mul(c);
    }
    b.print();
    cp_algo::checkpoint<1>();
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    solve();
}
