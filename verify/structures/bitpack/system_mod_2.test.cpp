// @brief System of Linear Equations (Mod 2)
#define PROBLEM "https://judge.yosupo.jp/problem/system_of_linear_equations_mod_2"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#include "cp-algo/structures/bit_array_util.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::structures;

const int maxn = (1 << 12) + 1;
bit_array<maxn> a[maxn];

void solve() {
    size_t n, m;
    cin >> n >> m;
    cp_algo::big_vector<cp_algo::big_string> As(n);
    for(size_t i = 0; i < n; i++) {
        cin >> As[i];
    }
    cp_algo::big_string bs;
    cin >> bs;
    for(size_t i = 0; i < n; i++) {
        As[i] += bs[i];
        from_string(a[i], As[i]);
    }
    cp_algo::big_vector<size_t> lead(n);
    auto vars = views::iota((size_t)0, m + 1);
    cp_algo::big_set<size_t> free(begin(vars), end(vars));
    for(size_t i = 0; i < n; i++) {
        for(size_t j = 0; j < i; j++) {
            if(a[i][lead[j]]) {
                a[i].xor_hint(a[j], lead[j]);
            }
        }
        lead[i] = ctz(a[i]);
        if(lead[i] == m) {
            cout << -1 << "\n";
            return;
        }
        if(lead[i] > m) {
            continue;
        }
        free.erase(lead[i]);
        for(size_t j = 0; j < i; j++) {
            if(a[j][lead[i]]) {
                a[j].xor_hint(a[i], lead[i]);
            }
        }
    }
    bit_array<maxn> x[maxn];
    for(auto [j, pj]: views::enumerate(free)) {
        x[j].set(pj);
        for(size_t i = 0; i < n; i++) {
            if(lead[i] < m && a[i][pj]) {
                x[j].set(lead[i]);
            }
        }
    }
    size_t rk = size(free) - 1;
    swap(x[0], x[rk]);
    cout << rk << "\n";
    for(size_t i = 0; i <= rk; i++) {
        cout << to_string(x[i]).substr(0, m) << "\n";
    }
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
