// @brief Matching on General Graph
#define PROBLEM "https://judge.yosupo.jp/problem/general_matching"
#include <bits/stdc++.h>
#pragma GCC optimize("Ofast,unroll-loops")
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/linalg/matrix.hpp"

using namespace std;
using namespace cp_algo::math;
using namespace cp_algo::linalg;
using namespace cp_algo::random;

const int64_t mod = 998244353;
using base = modint<mod>;

void solve() {
    int n, m;
    cin >> n >> m;
    matrix<base> T(n);
    for(int i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        base x = rng();
        T[u][v] += x;
        T[v][u] -= x;
    }
    auto [pivots, free] = matrix(T).echelonize();
    matrix<base> B(size(pivots));
    for(int i = 0; i < ssize(pivots); i++) {
        for(int j = 0; j < ssize(pivots); j++) {
            B[i][j] = T[pivots[i]][pivots[j]];
        }
    }
    B = B.inv().second;
    vector<pair<int, int>> ans;
    for(size_t i = 0; i < size(pivots); i++) {
        for(size_t j = 0; j < size(pivots); j++) {
            if(T[pivots[i]][pivots[j]] != 0 && B[i][j] != 0) {
                ans.emplace_back(pivots[i], pivots[j]);
                B.eliminate<gauss_mode::reverse>(i, j);
                B.eliminate<gauss_mode::reverse>(j, i);
                B.normalize();
                for(int k = 0; k < ssize(pivots); k++) {
                    B[i][k] = B[j][k] = 0;
                }
            }
        }
    }
    cout << size(ans) << "\n";
    for(auto [u, v]: ans) {
        cout << u << ' ' << v << "\n";
    }
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    // cin >> t;
    while(t--) {
        solve();
    }
}
