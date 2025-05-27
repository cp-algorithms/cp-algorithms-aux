// @brief Matrix Product (Mod 2)
#define PROBLEM "https://judge.yosupo.jp/problem/matrix_product_mod_2"
#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include "cp-algo/structures/bitpack.hpp"
#include "cp-algo/util/checkpoint.hpp"
#include <bits/stdc++.h>

using namespace std;
using cp_algo::structures::bitpack;

const int maxn = 1 << 12;
bitpack<maxn> a[maxn], b[maxn], c[maxn];
const int K = 8;
bitpack<maxn> precalc[1 << K];

void process_precalc(int i) {
    for(auto &it: precalc) {
        it = bitpack<maxn>();
    }
    for(int j = 0; j < K; j++) {
        int step = 1 << j;
        for(int k = 0; k < step; k++) {
            precalc[k + step] = precalc[k] ^ b[K * i + j];
        }
    }
}

void solve() {
    cp_algo::checkpoint("init");
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
    cp_algo::checkpoint("read");
    for(int j = 0; j < m; j += 64) {
        for(int z = 0; z < 64 / K; z++) {
            process_precalc(j / K + z);
            for(int i = 0; i < n; i++) {
                c[i] ^= precalc[uint8_t(a[i].word(j / 64) >> K * z)];
            }
        }
    }
    cp_algo::checkpoint("mul");
    for(int i = 0; i < n; i++) {
        row = c[i].to_string().substr(0, k);
        cout << row << "\n";
    }
    cp_algo::checkpoint("write");
    cp_algo::checkpoint<1>();
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
