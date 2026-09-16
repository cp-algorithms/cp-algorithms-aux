// @brief Number of Substrings
#define PROBLEM "https://judge.yosupo.jp/problem/number_of_substrings"
#include "cp-algo/structures/suffix_automaton.hpp"
#include <iostream>
#include "blazingio/blazingio.min.hpp"
using namespace cp_algo::structures;

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::string s;
    std::cin >> s;
    std::cout << suffix_automaton(std::move(s)).count_distinct() << '\n';
}
