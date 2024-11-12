
// @brief Matrix Product (Mod 2)
#define PROBLEM "https://judge.yosupo.jp/problem/matrix_product_mod_2"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("tune=native")
#include "cp-algo/structures/bitpack.hpp"
#include <bits/stdc++.h>

using namespace std;
using cp_algo::structures::bitpack;

const int maxn = 1 << 12;
bitpack<maxn> a[maxn], b[maxn], c[maxn];

void solve() {
    int n, m, k;
    cin >> n >> m >> k;
    string row;
    for(int i = 0; i < n; i++) {
        cin >> row;
        a[i] = row;
    }
    for(int i = 0; i < m; i++) {
        cin >> row;
        b[i] = row;
    }
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            if(a[i][j]) {
                c[i] ^= b[j];
            }
        }
    }
    for(int i = 0; i < n; i++) {
        row = c[i].to_string().substr(0, k);
        cout << row << "\n";
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
