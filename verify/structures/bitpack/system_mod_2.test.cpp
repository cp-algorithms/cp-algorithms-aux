// @brief System of Linear Equations (Mod 2)
// competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/system_of_linear_equations_mod_2
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/structures/bitpack.hpp"
#include <bits/stdc++.h>

using namespace std;
using cp_algo::structures::bitpack;

const int maxn = (1 << 12) + 1;
bitpack<maxn> a[maxn];

void solve() {
    size_t n, m;
    cin >> n >> m;
    vector<string> As(n);
    for(size_t i = 0; i < n; i++) {
        cin >> As[i];
    }
    string bs;
    cin >> bs;
    for(size_t i = 0; i < n; i++) {
        As[i] += bs[i];
        a[i] = As[i];
    }
    vector<size_t> lead(n);
    auto vars = views::iota((size_t)0, m + 1);
    set<size_t> free(begin(vars), end(vars));
    for(size_t i = 0; i < n; i++) {
        for(size_t j = 0; j < i; j++) {
            if(a[i][lead[j]]) {
                a[i].xor_hint(a[j], lead[j]);
            }
        }
        lead[i] = a[i].ctz();
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
    bitpack<maxn> x[maxn];
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
        cout << x[i].to_string().substr(0, m) << "\n";
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
