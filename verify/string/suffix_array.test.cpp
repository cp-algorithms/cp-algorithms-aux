// @brief Suffix Array
#define PROBLEM "https://judge.yosupo.jp/problem/suffixarray"
#pragma GCC target("popcnt")
#pragma GCC optimize("O3,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include "cp-algo/util/checkpoint.hpp"
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <numeric>
#include <vector>

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::string s;
    std::cin >> s;
    std::ranges::reverse(s);
    int n = (int)s.size();
    // Ordinary state i represents the prefix of length i; clones follow state n.
    struct clone_data { int len, pos; };
    std::vector<clone_data> clones;
    clones.reserve(n);
    cp_algo::big_vector<int> link(2 * n + 1);
    auto length = [&](int v) { return v <= n ? v : clones[v - n - 1].len; };
    auto position = [&](int v) { return v <= n ? v : clones[v - n - 1].pos; };
    // Pack two (letter + 1, state) pairs inline; bit 31 marks a dense row.
    // State IDs fit in 24 bits under the problem's n <= 500000 bound.
    struct transitions { uint32_t a = 0, b = 0; };
    cp_algo::big_vector<std::array<int, 26>> dense;
    // At most one dense row per state; avoid growth on branching structured strings.
    dense.reserve(link.size());
    cp_algo::big_vector<transitions> to;
    to.reserve(link.size());
    to.resize(n + 1);
    auto get = [&](int p, int x) {
        auto t = to[p];
        if(t.a & 0x80000000) {return dense[t.a & 0x7FFFFFFF][x];}
        if((t.a >> 24) == unsigned(x + 1)) {return int(t.a & 0xFFFFFF);}
        if((t.b >> 24) == unsigned(x + 1)) {return int(t.b & 0xFFFFFF);}
        return 0;
    };
    auto set = [&](int p, int x, int v) {
        auto &t = to[p];
        auto encoded = uint32_t((x + 1) << 24) | v;
        if(t.a & 0x80000000) {dense[t.a & 0x7FFFFFFF][x] = v;}
        else if(!t.a || (t.a >> 24) == unsigned(x + 1)) {t.a = encoded;}
        else if(!t.b || (t.b >> 24) == unsigned(x + 1)) {t.b = encoded;}
        else {
            std::array<int, 26> row{};
            row[(t.a >> 24) - 1] = t.a & 0xFFFFFF;
            row[(t.b >> 24) - 1] = t.b & 0xFFFFFF;
            row[x] = v;
            t.a = 0x80000000 | uint32_t(dense.size());
            dense.push_back(row);
        }
    };
    cp_algo::checkpoint("init");
    int last = 0;
    for(char c: s) {
        int x = c - 'a', p = last;
        ++last;
        for(; !get(p, x); p = link[p]) {
            set(p, x, last);
        }
        int q = get(p, x);
        if(q != last) {
            if(length(q) == length(p) + 1) {
                link[last] = q;
            } else {
                int clone = n + 1 + (int)clones.size();
                clones.push_back({length(p) + 1, position(q)});
                link[clone] = link[q];
                auto row = to[q];
                if(row.a & 0x80000000) {
                    auto copy = dense[row.a & 0x7FFFFFFF];
                    row.a = 0x80000000 | uint32_t(dense.size());
                    dense.push_back(copy);
                }
                to.push_back(row);
                link[last] = link[q] = clone;
                for(; get(p, x) == q; p = link[p]) {
                    set(p, x, clone);
                }
            }
        }
    }
    int size = n + 1 + (int)clones.size();
    dense = decltype(dense){};
    cp_algo::checkpoint("build");
    // Children have distinct incoming letters. Reuse b for their letter masks.
    cp_algo::big_vector<int> offset(size + 1);
    for(auto &t: to) {
        t.b = 0;
    }
    for(int i = 1; i < size; i++) {
        int p = link[i], x = s[position(i) - length(p) - 1] - 'a';
        to[p].b |= 1U << x;
        link[i] = (x << 24) | p;
    }
    for(int i = 0; i < size; i++) {
        offset[i + 1] = offset[i] + __builtin_popcount(to[i].b);
    }
    // Count smaller letters to scatter child IDs into a in lexicographic order.
    for(int i = 1; i < size; i++) {
        int p = link[i] & 0xFFFFFF, x = unsigned(link[i]) >> 24;
        int rank = __builtin_popcount(to[p].b & ((1U << x) - 1));
        to[offset[p] + rank].a = i;
    }
    cp_algo::checkpoint("tree");
    // The masks and suffix links are now dead; reuse them for stack/output.
    int top = 1, count = 0;
    to[0].b = 0;
    while(top) {
        int u = to[--top].b;
        if(u && u <= n) {
            link[count++] = n - u;
        }
        for(int j = offset[u + 1]; j > offset[u];) {
            to[top++].b = to[--j].a;
        }
    }
    cp_algo::checkpoint("dfs");
    for(int i = 0; i < n; i++) {
        std::cout << link[i] << ' ';
    }
    std::cout << '\n';
    cp_algo::checkpoint("write");
    cp_algo::checkpoint<1>();
}
