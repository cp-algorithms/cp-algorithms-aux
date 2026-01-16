// @brief Chromatic Number
#define PROBLEM "https://judge.yosupo.jp/problem/chromatic_number"
#pragma GCC optimize("O3,unroll-loops")
#include <bits/allocator.h>
#pragma GCC target("avx2")
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#define CP_ALGO_CHECKPOINT
#include "cp-algo/number_theory/modint.hpp"
#include "cp-algo/math/subset_convolution.hpp"
#include "cp-algo/math/poly.hpp"
#include <bits/stdc++.h>

using namespace std;

const int mod = 998244353;
using base = cp_algo::math::modint<mod>;
using polyn = cp_algo::math::poly_t<base>;

cp_algo::big_vector<base> indep_subsets(auto const& adj) {
    uint32_t n = (uint32_t)size(adj);
    uint32_t masks = 1 << n;
    cp_algo::big_vector<base> indep(masks);
    indep[0] = base(1);
    for(uint32_t v = 0; v < n; v++) {
        uint32_t delta = 1 << v;
        if (adj[v] & delta) continue;
        for(uint32_t mask = 0; mask < delta; mask++) {
            indep[mask + delta] = (adj[v] & mask) ? base(0) : indep[mask];
        }
    }
    return indep;
}

void solve() {
    size_t n, m;
    cin >> n >> m;
    vector<uint64_t> adj(n);
    for(size_t i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        adj[u] |= 1 << v;
        adj[v] |= 1 << u;
    }
    auto indep = indep_subsets(adj);
    cp_algo::big_vector<base> w(1 << n);
    w.back() = 1; // w[S] = 1 if S is the full set, else 0
    auto Y = cp_algo::math::subset_power_projection<base>(indep, w, n+1);
    size_t ans = 0;
    while(Y[ans] == base(0)) {
        ans++;
    }
    cout << ans << "\n";
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
