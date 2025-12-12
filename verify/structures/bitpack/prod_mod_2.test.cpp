// @brief Matrix Product (Mod 2)
#define PROBLEM "https://judge.yosupo.jp/problem/matrix_product_mod_2"
#pragma GCC optimize("O3,unroll-loops")
#define CP_ALGO_CHECKPOINT
#pragma GCC target("avx2")
#include "cp-algo/structures/bit_array_util.hpp"
#include "cp-algo/util/checkpoint.hpp"
#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo::structures;

const int maxn = 1 << 12;
const size_t K = 8;

using bitpack = cp_algo::structures::bit_array<maxn>;

bitpack a[maxn], b[maxn], c[maxn];
bitpack precalc[1 << K];

void process_precalc(int i) {
    for(size_t j = 0; j < K; j++) {
        int step = 1 << j;
        for(int k = 0; k < step; k++) {
            precalc[k + step] = precalc[k] ^ b[i + j];
        }
    }
}

void solve() {
    int n, m, k;
    cin >> n >> m >> k;
    cp_algo::checkpoint("init");
    string row;
    for(int i = 0; i < n; i++) {
        cin >> row;
        from_string(a[i], row);
    }
    for(int i = 0; i < m; i++) {
        cin >> row;
        from_string(b[i], row);
    }
    for(auto &it: c) {
        it.resize(k);
    }
    for(auto &it: precalc) {
        it.resize(k);
    }
    cp_algo::checkpoint("read");
    const int width = bitpack::width;
    for(int j = 0; j < m; j += width) {
        for(int offset = 0; offset < width; offset += K) {
            process_precalc(j + offset);
            for(int i = 0; i < n; i++) {
                c[i] ^= precalc[uint8_t(a[i].word(j / width) >> offset)];
            }
        }
    }
    cp_algo::checkpoint("mul");
    for(int i = 0; i < n; i++) {
        cout << to_string(c[i]) << "\n";
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
