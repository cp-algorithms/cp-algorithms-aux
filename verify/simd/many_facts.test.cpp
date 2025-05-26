// @brief Many Factorials
#define PROBLEM "https://judge.yosupo.jp/problem/many_factorials"
#pragma GCC optimize("Ofast,unroll-loops")
#include <bits/stdc++.h>
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/simd.hpp"
#include "cp-algo/math/common.hpp"

using namespace std;
using namespace cp_algo;

constexpr int mod = 998244353;
constexpr auto mod4 = u64x4() + mod;
constexpr auto imod4 = u64x4() - math::inv2(mod);

void facts_inplace(vector<int> &args) {
    constexpr int block = 1 << 16;
    static basic_string<size_t> args_per_block[mod / block];
    uint64_t limit = 0;
    for(auto [i, x]: args | views::enumerate) {
        if(x < mod / 2) {
            limit = max(limit, uint64_t(x));
            args_per_block[x / block].push_back(i);
        } else {
            limit = max(limit, uint64_t(mod - x - 1));
            args_per_block[(mod - x - 1) / block].push_back(i);
        }
    }
    uint64_t b2x32 = (1ULL << 32) % mod;
    uint64_t fact = 1;
    for(uint64_t b = 0; b <= limit; b += 4 * block) {
        u64x4 cur[4];
        static array<u64x4, block / 4> prods[4];
        for(int z = 0; z < 4; z++) {
            for(int j = 0; j < 4; j++) {
                cur[z][j] = b + z * block + j * block / 4;
                prods[z][0][j] = cur[z][j] + !(b || z || j);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
                cur[z][j] = cur[z][j] * b2x32 % mod;
#pragma GCC diagnostic pop
            }
        }
        for(int i = 1; i < block / 4; i++) {
            for(int z = 0; z < 4; z++) {
                cur[z] += b2x32;
                cur[z] = cur[z] >= mod ? cur[z] - mod : cur[z];
                prods[z][i] = montgomery_mul(prods[z][i - 1], cur[z], mod4, imod4);
            }
        }
        for(int z = 0; z < 4; z++) {
            uint64_t bl = b + z * block;
            for(auto i: args_per_block[bl / block]) {
                size_t x = args[i];
                if(x >= mod / 2) {
                    x = mod - x - 1;
                }
                x -= bl;
                auto pre_blocks = x / (block / 4);
                auto in_block = x % (block / 4);
                auto ans = fact * prods[z][in_block][pre_blocks] % mod;
                for(size_t j = 0; j < pre_blocks; j++) {
                    ans = ans * prods[z].back()[j] % mod;
                }
                if(args[i] >= mod / 2) {
                    ans = math::bpow(ans, mod - 2, 1ULL, [](auto a, auto b){return a * b % mod;});
                    args[i] = int(x % 2 ? ans : mod - ans);
                } else {
                    args[i] = int(ans);
                }
            }
            args_per_block[bl / block].clear();
            for(int j = 0; j < 4; j++) {
                fact = fact * prods[z].back()[j] % mod;
            }
        }
    }
}

void solve() {
    int n;
    cin >> n;
    vector<int> args(n);
    for(auto &x : args) {cin >> x;}
    facts_inplace(args);
    for(auto it: args) {cout << it << "\n";}
}

signed main() {
    //freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    int t = 1;
    //cin >> t;
    while(t--) {
        solve();
    }
}
