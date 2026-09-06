// @brief Shift of Sampling Points of Polynomial
#define PROBLEM "https://judge.yosupo.jp/problem/shift_of_sampling_points_of_polynomial"
#define CP_ALGO_MAXN 1 << 20
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/combinatorics.hpp"
#include "cp-algo/math/poly/base.hpp"

using namespace std;
using namespace cp_algo::math;

const int mod = 998244353;
using base = modint<mod>;
using polyn = poly_t<base>;

void solve() {
    int n, m, c;
    cin >> n >> m >> c;
    polyn::Vector a(n);
    for(auto &it: a) {cin >> it;}

    // f(x) = prod_j(x-j) * sum_i (-1)^(n-1-i) f(i) / (i! (n-1-i)! (x-i)).
    // The sums for consecutive x are one convolution with coefficients 1/j of log(1/(1-x)).
    auto weighted = polyn::Vector(std::from_range, views::iota(0, n) | views::transform([&](int i) {
        base v = a[i] * rfact<base>(i) * rfact<base>(n - 1 - i);
        return (n - 1 - i) & 1 ? -v : v;
    }));
    auto denominators = views::iota(c - n + 1, c + m) | views::transform([](int x) {
        return base(x) == base(0) ? base(1) : base(x);
    });
    auto log_segment = bulk_invs<base>(denominators);
    for(int z: {n - 1 - c, n - 1 - c + mod}) {
        if(0 <= z && z < ssize(log_segment)) {log_segment[z] = 0;}
    }
    auto shifted = polyn(std::move(weighted)) * polyn(log_segment);

    base product = ranges::fold_left(views::iota(0, n), base(1), [c](base p, int i) {return p * base(c - i);});
    for(int i = 0, x = c; i < m; i++, x = (x + 1) % mod) {
        cout << (x < n ? a[x] : shifted[n - 1 + i] * product) << ' ';
        base next = base(x) + 1;
        product = next == base(n) ? fact<base>(n) : product * next * log_segment[i];
    }
    cout << '\n';
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t;
    t = 1;// cin >> t;
    while(t--) {
        solve();
    }
}
