
// @brief Determinant of Matrix (Mod 2)
#define PROBLEM "https://judge.yosupo.jp/problem/matrix_det_mod_2"
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("avx2,tune=native")
#include "cp-algo/data_structures/bitpack.hpp"
#include <bits/stdc++.h>

using namespace std;
using cp_algo::data_structures::bitpack;

const int maxn = 1 << 12;
bitpack<maxn> a[maxn];

void solve() {
    int n;
    cin >> n;
    string row;
    vector<int> lead(n);
    for(int i = 0; i < n; i++) {
        cin >> row;
        a[i] = row;
        for(int j = 0; j < i; j++) {
            if(a[i][lead[j]]) {
                a[i] ^= a[j];
            }
        }
        lead[i] = a[i].ctz();
        if(lead[i] == maxn) {
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
