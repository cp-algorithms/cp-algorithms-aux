// @brief Static Convex Hull
#define PROBLEM "https://judge.yosupo.jp/problem/static_convex_hull"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/geometry/convex_hull.hpp"
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
    auto res = convex_hull(r);
    cout << size(res) << "\n";
    for(auto it: res) {
        it.print();
    }
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
