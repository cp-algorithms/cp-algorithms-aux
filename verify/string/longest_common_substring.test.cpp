// @brief Longest Common Substring
#define PROBLEM "https://judge.yosupo.jp/problem/longest_common_substring"
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#define CP_ALGO_CHECKPOINT
#include "cp-algo/util/checkpoint.hpp"
#include <bits/stdc++.h>
#include "cp-algo/util/big_alloc.hpp"
#include <iostream>
#include "blazingio/blazingio.min.hpp"

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::string s, t;
    std::cin >> s >> t;
    bool swapped = s.size() > t.size();
    if(swapped) std::swap(s, t);
    std::array<int, 26> code;
    code.fill(-1);
    int alphabet = 0;
    for(auto str: {&s, &t}) {
        for(char &c: *str) {
            int &rank = code[c - 'a'];
            if(rank < 0) rank = alphabet++;
            c = char(rank);
        }
    }
    cp_algo::big_vector<int> len(2 * s.size() + 1), link(len.size()), pos(len.size());
    cp_algo::big_vector<int> transitions(len.size() * alphabet);
    auto edge = [&](int v, int c) -> int& { return transitions[size_t(v) * alphabet + c]; };
    cp_algo::checkpoint("init");
    int last = 0, size = 1;
    for(char c: s) {
        int x = c, p = last;
        last = size++;
        len[last] = pos[last] = len[p] + 1;
        for(; !edge(p, x); p = link[p]) {
            edge(p, x) = last;
        }
        int q = edge(p, x);
        if(q != last) {
            if(len[q] == len[p] + 1) {
                link[last] = q;
            } else {
                int clone = size++;
                link[clone] = link[q];
                pos[clone] = pos[q];
                std::copy_n(&edge(q, 0), alphabet, &edge(clone, 0));
                len[clone] = len[p] + 1;
                link[last] = link[q] = clone;
                for(; edge(p, x) == q; p = link[p]) {
                    edge(p, x) = clone;
                }
            }
        }
    }
    cp_algo::checkpoint("build");
    int v = 0, length = 0, best = 0, end_s = 0, end_t = 0;
    for(int i = 0; i < (int)t.size(); i++) {
        int x = t[i];
        while(v && !edge(v, x)) {
            v = link[v];
            length = len[v];
        }
        v = edge(v, x);
        length = v ? length + 1 : 0;
        if(length > best) {
            best = length;
            end_s = pos[v];
            end_t = i + 1;
        }
    }
    cp_algo::checkpoint("query");
    cp_algo::checkpoint<1>();
    if(swapped) std::swap(end_s, end_t);
    std::cout << end_s - best << ' ' << end_s << ' '
              << end_t - best << ' ' << end_t << '\n';
}
