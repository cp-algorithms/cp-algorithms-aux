// @brief Longest Common Substring
#define PROBLEM "https://judge.yosupo.jp/problem/longest_common_substring"
#include <bits/stdc++.h>
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#define CP_ALGO_CHECKPOINT
#include "cp-algo/structures/suffix_automaton.hpp"
#include <iostream>
#include "blazingio/blazingio.min.hpp"
using namespace cp_algo::structures;

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::string s, t;
    std::cin >> s >> t;
    bool swapped = s.size() > t.size();
    if(swapped) std::swap(s, t);
    auto [a, b, c, d] = suffix_automaton(std::move(s)).longest_common_substring(t);
    if(swapped) { std::swap(a, c); std::swap(b, d); }
    std::cout << a << ' ' << b << ' ' << c << ' ' << d << '\n';
    cp_algo::checkpoint<1>();
}
