// @brief Many Factorials
#define PROBLEM "https://judge.yosupo.jp/problem/many_factorials"
#pragma GCC optimize("Ofast,unroll-loops")
#include <bits/stdc++.h>
//#define CP_ALGO_CHECKPOINT
#include "blazingio/blazingio.min.hpp"
#include "cp-algo/util/checkpoint.hpp"
#include "cp-algo/util/simd.hpp"
#include "cp-algo/util/bump_alloc.hpp"
#include "cp-algo/math/common.hpp"

using namespace std;
using namespace cp_algo;

constexpr int mod = 998244353;
constexpr int imod = -math::inv2(mod);

template<int maxn = 100'000>
vector<int> facts(vector<int> const& args) {
    constexpr int accum = 4;
    constexpr int simd_size = 8;
    constexpr int block = 1 << 18;
    constexpr int subblock = block / simd_size;
    using T = array<int, 2>;
    using alloc = bump_alloc<T, 30 * maxn>;
    basic_string<T, char_traits<T>, alloc> odd_args_per_block[mod / subblock];
    basic_string<T, char_traits<T>, alloc> reg_args_per_block[mod / subblock];
    constexpr int limit_reg = mod / 64;
    int limit_odd = 0;

    vector<int> res(size(args), 1);
    auto prod_mod = [&](uint64_t a, uint64_t b) {
        return (a * b) % mod;
    };
    for(auto [i, xy]: views::zip(args, res) | views::enumerate) {
        auto [x, y] = xy;
        auto t = x;
        if(t >= mod / 2) {
            t = mod - t - 1;
            y = t % 2 ? 1 : mod - 1;
        }
        int pw = 0;
        while(t > limit_reg) {
            limit_odd = max(limit_odd, (t - 1) / 2);
            odd_args_per_block[(t - 1) / 2 / subblock].push_back({int(i), (t - 1) / 2});
            t /= 2;
            pw += t;
        }
        reg_args_per_block[t / subblock].push_back({int(i), t});
        y = int(y * math::bpow(2, pw, 1ULL, prod_mod) % mod);
    }
    cp_algo::checkpoint("init");
    uint32_t b2x32 = (1ULL << 32) % mod;
    auto process = [&](int limit, auto &args_per_block, auto step, auto &&proj) {
        uint64_t fact = 1;
        for(int b = 0; b <= limit; b += accum * block) {
            u32x8 cur[accum];
            static array<u32x8, subblock> prods[accum];
            for(int z = 0; z < accum; z++) {
                for(int j = 0; j < simd_size; j++) {
                    cur[z][j] = uint32_t(b + z * block + j * subblock);
                    cur[z][j] = proj(cur[z][j]);
                    prods[z][0][j] = cur[z][j] + !cur[z][j];
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
                    cur[z][j] = uint32_t(uint64_t(cur[z][j]) * b2x32 % mod);
    #pragma GCC diagnostic pop
                }
            }
            for(int i = 1; i < block / simd_size; i++) {
                for(int z = 0; z < accum; z++) {
                    cur[z] += step;
                    cur[z] = cur[z] >= mod ? cur[z] - mod : cur[z];
                    prods[z][i] = montgomery_mul(prods[z][i - 1], cur[z], mod, imod);
                }
            }
            cp_algo::checkpoint("inner loop");
            for(int z = 0; z < accum; z++) {
                for(int j = 0; j < simd_size; j++) {
                    int bl = b + z * block + j * subblock;
                    for(auto [i, x]: args_per_block[bl / subblock]) {
                        auto ans = fact * prods[z][x - bl][j] % mod;
                        res[i] = int(res[i] * ans % mod);
                    }
                    fact = fact * prods[z].back()[j] % mod;
                }
            }
            cp_algo::checkpoint("mul ans");
        }
    };
    uint32_t b2x33 = (1ULL << 33) % mod;
    process(limit_reg, reg_args_per_block, b2x32, identity{});
    process(limit_odd, odd_args_per_block, b2x33, [](uint32_t x) {return 2 * x + 1;});
    for(auto [i, x]: res | views::enumerate) {
        if (args[i] >= mod / 2) {
            x = int(math::bpow(x, mod - 2, 1ULL, prod_mod));
        }
    }
    cp_algo::checkpoint("inv ans");
    return res;
}

void solve() {
    int n;
    cin >> n;
    vector<int> args(n);
    for(auto &x : args) {cin >> x;}
    cp_algo::checkpoint("read");
    auto res = facts(args);
    for(auto it: res) {cout << it << "\n";}
    cp_algo::checkpoint("write");
    cp_algo::checkpoint<1>();
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
