// @brief Multidimensional Convolution (Truncated)
#define PROBLEM "https://judge.yosupo.jp/problem/multivariate_convolution"
#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#define CP_ALGO_CHECKPOINT
#pragma GCC target("avx2")
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/io.hpp"
#include "cp-algo/math/multivar.hpp"

using cp_algo::io::read_reduced;
using cp_algo::io::read;

using namespace std;
using namespace cp_algo::math::fft;

const int mod = 998244353;
using base = cp_algo::math::modint<mod>;

void solve() {
    int k;
    cin >> k;
    cp_algo::big_vector<size_t> N(k);
    read(N);
    multivar<base> a(N), b(N);
    read_reduced(a.data);
    cp_algo::checkpoint("multivar read");
    read_reduced(b.data);
    cp_algo::checkpoint("multivar read");
    a.mul(b);
    a.print();
    cp_algo::checkpoint<1>();
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    solve();
}
