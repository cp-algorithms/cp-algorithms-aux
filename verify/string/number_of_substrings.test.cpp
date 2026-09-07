// @brief Number of Substrings
#define PROBLEM "https://judge.yosupo.jp/problem/number_of_substrings"
#include <iostream>
#include "blazingio/blazingio.min.hpp"
#include <bits/stdc++.h>
#include <sys/mman.h>

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::string s;
    std::cin >> s;
    struct node {
        int len, link;
        std::array<int, 26> to;
    };
    // Zero-filled pages are populated as states are used; retain huge-page support.
    size_t bytes = (2 * s.size() + 1) * sizeof(node);
    auto nodes = static_cast<node*>(mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                                        MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
    assert(nodes != MAP_FAILED);
    madvise(nodes, bytes, MADV_HUGEPAGE);
    int last = 0, size = 1;
    int64_t answer = 0;
    for(char c: s) {
        int x = c - 'a', p = last;
        last = size++;
        nodes[last].len = nodes[p].len + 1;
        for(; !nodes[p].to[x]; p = nodes[p].link) {
            nodes[p].to[x] = last;
        }
        int q = nodes[p].to[x];
        if(q != last) {
            if(nodes[q].len == nodes[p].len + 1) {
                nodes[last].link = q;
            } else {
                int clone = size++;
                nodes[clone] = nodes[q];
                nodes[clone].len = nodes[p].len + 1;
                nodes[last].link = nodes[q].link = clone;
                for(; nodes[p].to[x] == q; p = nodes[p].link) {
                    nodes[p].to[x] = clone;
                }
            }
        }
        answer += nodes[last].len - nodes[nodes[last].link].len;
    }
    std::cout << answer << '\n';
}
