#include <cassert>
#include <bits/stdc++.h>
#ifdef TEST_MINIFIED
#include "cp-algo/min/structures/eertree.hpp"
#else
#include "cp-algo/structures/eertree.hpp"
#endif
template <int Sigma, char First> void check(std::string const &text) {
    cp_algo::structures::eertree<Sigma, First> tree(text.size());
    std::map<std::string, int> ids;
    for(size_t end = 0; end < text.size(); ++end) {
        tree.add_letter(text[end]);
        for(size_t start = 0; start <= end; ++start) {
            auto word = text.substr(start, end - start + 1);
            if(std::equal(word.begin(), word.end(), word.rbegin())) {
                if(!ids.contains(word)) {
                    int id = ids.size() + 2;
                    ids[word] = id;
                }
                assert(tree.sufpal() == ids.at(word));
                break;
            }
        }
    }
    std::ostringstream output;
    auto previous = std::cout.rdbuf(output.rdbuf());
    tree.print();
    std::cout.rdbuf(previous);
    std::istringstream input(output.str());
    size_t count;
    input >> count;
    assert(count == ids.size());
}
int main() {
    check<10, '0'>("000000");
    check<2, char(127)>(std::string{char(127), char(128), char(127), char(128)});
    std::mt19937 rng(2931);
    for(int test = 0; test < 500; ++test) {
        std::string digits(1 + rng() % 80, '0');
        std::string upper(digits.size(), 'A'), lower(digits.size(), 'a');
        for(size_t i = 0; i < digits.size(); ++i) {
            digits[i] += rng() % 10;
            upper[i] += rng() % 26;
            lower[i] += rng() % 26;
        }
        check<10, '0'>(digits);
        check<26, 'A'>(upper);
        check<26, 'a'>(lower);
        std::string bytes(digits.size(), 0);
        for(char &c : bytes)
            c = char(rng() % 255);
        check<255, 0>(bytes);
    }
}
