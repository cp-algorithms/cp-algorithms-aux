// @brief Closest Pair of Points
#define PROBLEM "https://judge.yosupo.jp/problem/closest_pair"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/geometry/closest_pair.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::geometry;
using point = point_t<int64_t>;

void solve() {
    int n;
    cin >> n;
    cp_algo::big_vector<point> r(n);
    for(auto &it: r) {
        it.read();
    }
    auto [a, b] = closest_pair(r);
    cout << a << ' ' << b << "\n";
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    cin >> t;
    while(t--) {
        solve();
    }
}
