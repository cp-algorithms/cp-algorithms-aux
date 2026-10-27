// @brief Multipoint Evaluation
#define PROBLEM "https://judge.yosupo.jp/problem/multipoint_evaluation"
#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/io.hpp"
#include "cp-algo/math/poly/eval.hpp"

using cp_algo::io::read_reduced;

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n, m;
    cin >> n >> m;
    polyn::Vector f(n), x(m);
    read_reduced(f);
    read_reduced(x);
    polyn(eval(polyn(std::move(f)), x)).print(m);
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
