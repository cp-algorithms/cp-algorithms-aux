// @brief Suffix Array
#define PROBLEM "https://judge.yosupo.jp/problem/suffixarray"
#pragma GCC target("popcnt")
#pragma GCC optimize("O3,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include "cp-algo/structures/suffix_automaton.hpp"
#include <iostream>
#include "blazingio/blazingio.min.hpp"
using namespace cp_algo::structures;

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::string s;
    std::cin >> s;
    for(int i: suffix_array(std::move(s))) std::cout << i << ' ';
    std::cout << '\n';
    cp_algo::checkpoint("write");
    cp_algo::checkpoint<1>();
}
