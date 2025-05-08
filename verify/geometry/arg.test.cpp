// @brief Sort Points by Argument
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/sort_points_by_argument
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("tune=native")
#include "cp-algo/geometry/point.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::geometry;
using point = point_t<int64_t>;

void solve() {
    int n;
    cin >> n;
    vector<point> pts(n);
    for(auto &r: pts) {
        r.read();
    }
    ranges::sort(pts, point::ccw_abs);
    for(auto r: pts) {
        r.print();
    }
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
