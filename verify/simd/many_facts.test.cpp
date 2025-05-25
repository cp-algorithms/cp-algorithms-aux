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
    for(uint64_t b = 0; b <= limit; b += block) {
        u64x4 cur = {b, b + block / 4, b + block / 2, b + 3 * block / 4};
        static array<u64x4, block / 4> prods;
        prods[0] = u64x4{cur[0] + !b, cur[1], cur[2], cur[3]};
        cur = cur * b2x32 % mod;
        for(int i = 1; i < block / 4; i++) {
            cur += b2x32;
            cur = cur >= mod ? cur - mod : cur;
            prods[i] = montgomery_mul(prods[i - 1], cur, mod4, imod4);
        }
        for(auto i: args_per_block[b / block]) {
            size_t x = args[i];
            if(x >= mod / 2) {
                x = mod - x - 1;
            }
            x -= b;
            auto pre_blocks = x / (block / 4);
            auto in_block = x % (block / 4);
            auto ans = fact * prods[in_block][pre_blocks] % mod;
            for(size_t z = 0; z < pre_blocks; z++) {
                ans = ans * prods.back()[z] % mod;
            }
            if(args[i] >= mod / 2) {
                ans = math::bpow(ans, mod - 2, 1ULL, [](auto a, auto b){return a * b % mod;});
                args[i] = int(x % 2 ? ans : mod - ans);
            } else {
                args[i] = int(ans);
            }
        }
        args_per_block[b / block].clear();
        for(int z = 0; z < 4; z++) {
            fact = fact * prods.back()[z] % mod;
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
