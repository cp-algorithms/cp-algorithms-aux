
// @brief Determinant of Matrix (Mod 2)
#define PROBLEM "https://judge.yosupo.jp/problem/matrix_det_mod_2"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/structures/bitpack.hpp"
#include <bits/stdc++.h>

using namespace std;
using cp_algo::structures::bitpack;

const int maxn = 1 << 12;
bitpack<maxn> a[maxn];

void solve() {
    size_t n;
    cin >> n;
    string row;
    vector<size_t> lead(n);
    for(size_t i = 0; i < n; i++) {
        cin >> row;
        a[i] = row;
        for(size_t j = 0; j < i; j++) {
            if(a[i][lead[j]]) {
                a[i].xor_hint(a[j], lead[j]);
            }
        }
        lead[i] = a[i].ctz();
        if(lead[i] == n) {
            cout << 0 << "\n";
            return;
        }
    }
    cout << 1 << "\n";
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
