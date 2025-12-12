// @brief Counting Eulerian Circuits
#define PROBLEM "https://judge.yosupo.jp/problem/counting_eulerian_circuits"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include <bits/stdc++.h>
//#include "blazingio/blazingio.min.hpp"
#include "cp-algo/math/combinatorics.hpp"
#include "cp-algo/linalg/matrix.hpp"

using namespace std;
using namespace cp_algo::math;
using namespace cp_algo::linalg;

const int64_t mod = 998244353;
using base = modint<mod>;

void solve() {
    int n, m;
    cin >> n >> m;
    matrix<base> a(n);
    int r = 0;
    cp_algo::big_vector<int> indeg(n), outdeg(n);
    for(int i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        a[u][v] -= 1;
        a[v][v] += 1;
        outdeg[u]++;
        indeg[v]++;
        r = v;
    }
    a[r][r] = fact<base>(indeg[r] - 1);
    for(int i = 0; i < n; i++) {
        if(i == r) {
            continue;
        }
        if(indeg[i] != outdeg[i]) {
            cout << 0 << "\n";
            return;
        }
        if(indeg[i] != 0) {
            a[i] *= fact<base>(indeg[i] - 1);
        } else {
            a[i][i] = 1;
        }
        a[r][i] = a[i][r] = 0;
    }
    cout << a.det() << "\n";
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
