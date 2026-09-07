// @brief Suffix Array
#define PROBLEM "https://judge.yosupo.jp/problem/suffixarray"
#pragma GCC optimize("O3,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include "cp-algo/util/checkpoint.hpp"
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>
#include <vector>

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::string s;
    std::cin >> s;
    std::ranges::reverse(s);
    int n = (int)s.size();
    cp_algo::big_vector<int> len(2 * n + 1), link(2 * n + 1), pos(2 * n + 1);
    using transitions = std::array<int, 26>;
    cp_algo::big_vector<transitions> to(len.size());
    cp_algo::checkpoint("init");
    int last = 0, size = 1;
    for(char c: s) {
        int x = c - 'a', p = last;
        last = size++;
        len[last] = pos[last] = len[p] + 1;
        for(; !to[p][x]; p = link[p]) {
            to[p][x] = last;
        }
        int q = to[p][x];
        if(q != last) {
            if(len[q] == len[p] + 1) {
                link[last] = q;
            } else {
                int clone = size++;
                link[clone] = link[q];
                pos[clone] = pos[q];
                to[clone] = to[q];
                len[clone] = len[p] + 1;
                link[last] = link[q] = clone;
                for(; to[p][x] == q; p = link[p]) {
                    to[p][x] = clone;
                }
            }
        }
    }
    cp_algo::checkpoint("build");
    // Store the suffix-link tree compactly, with each node's children in letter order.
    struct edge { int child; char letter; };
    cp_algo::big_vector<int> offset(size + 1);
    cp_algo::big_vector<edge> edges(size - 1);
    for(int i = 1; i < size; i++) {
        offset[link[i]]++;
    }
    std::partial_sum(offset.begin(), offset.end(), offset.begin());
    for(int i = 1; i < size; i++) {
        int p = link[i];
        edges[--offset[p]] = {i, s[pos[i] - len[p] - 1]};
    }
    for(int i = 0; i < size; i++) {
        std::sort(edges.begin() + offset[i], edges.begin() + offset[i + 1],
                  [](edge a, edge b) { return a.letter < b.letter; });
    }
    cp_algo::checkpoint("tree");
    std::vector<int> stack{0}, answer;
    answer.reserve(n);
    while(!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        if(u && len[u] == pos[u]) {
            answer.push_back(n - pos[u]);
        }
        for(int j = offset[u + 1]; j > offset[u];) {
            stack.push_back(edges[--j].child);
        }
    }
    cp_algo::checkpoint("dfs");
    for(int i: answer) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
    cp_algo::checkpoint("write");
    cp_algo::checkpoint<1>();
}
