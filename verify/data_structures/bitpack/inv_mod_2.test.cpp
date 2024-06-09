// @brief Inverse Matrix (Mod 2)
#define PROBLEM "https://judge.yosupo.jp/problem/inverse_matrix_mod_2"
#pragma GCC optimize("Ofast,unroll-loops")
#include "cp-algo/data_structures/bitpack.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::data_structures;

const int maxn = 1 << 13;
bitpack<maxn> a[maxn];

void solve() {
    int n;
    cin >> n;
    string row;
    vector<int> lead(n);
    for(int i = 0; i < n; i++) {
        cin >> row;
        a[i] = row;
        a[i].set(n + i);
        for(int j = 0; j < i; j++) {
            if(a[i][lead[j]]) {
                a[i].xor_hint(a[j], lead[j]);
            }
        }
        lead[i] = a[i].ctz();
        if(lead[i] >= n) {
            cout << -1 << "\n";
            return;
        }
        for(int j = 0; j < i; j++) {
            if(a[j][lead[i]]) {
                a[j].xor_hint(a[i], lead[i]);
            }
        }
    }
    for(int i = 0; i < n; i++) {
        while(lead[i] != i) {
            swap(a[i], a[lead[i]]);
            swap(lead[i], lead[lead[i]]);
        }
        cout << a[i].to_string().substr(n, n) << "\n";
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
