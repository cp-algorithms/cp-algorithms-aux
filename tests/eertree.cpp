#include <bits/stdc++.h>
#include <cassert>
#include "cp-algo/structures/eertree.hpp"
using namespace cp_algo::structures;

bool palindrome(std::string const& s) {
    return std::equal(s.begin(), s.end(), s.rbegin());
}
void check(std::string const& input) {
    eertree tree(input.size());
    std::vector<std::string> word(2);
    std::set<std::string> distinct;
    for(size_t end = 1; end <= input.size(); end++) {
        tree.add_letter(input[end - 1]);
        std::string longest;
        for(size_t begin = 0; begin < end; begin++) {
            auto w = input.substr(begin, end - begin);
            if(palindrome(w)) {
                distinct.insert(w);
                if(w.size() > longest.size()) longest = w;
            }
        }
        int id = tree.sufpal();
        if(id == (int)word.size()) word.push_back(longest);
        assert(id >= 2 && id < (int)word.size() && word[id] == longest);
    }
    std::ostringstream out;
    for(int c: {-256, 256, 257, 511, 1000}) {
        assert(tree.get(0, c) == 0 && tree.get(1, c) == 0);
    }
    auto saved = std::cout.rdbuf(out.rdbuf());
    tree.print();
    std::cout.rdbuf(saved);
    std::istringstream in(out.str());
    int count;
    in >> count;
    assert(count == (int)distinct.size() && count + 2 == (int)word.size());
    for(int i = 2; i < count + 2; i++) {
        int parent, link;
        in >> parent >> link;
        if(word[i].size() == 1) assert(parent == 1);
        else assert(word[parent] == word[i].substr(1, word[i].size() - 2));
        std::string suffix;
        for(size_t begin = 1; begin < word[i].size(); begin++) {
            auto w = word[i].substr(begin);
            if(palindrome(w)) {suffix = w; break;}
        }
        assert(word[link] == suffix);
    }
}
int main() {
    check("");
    std::mt19937 rng(73159);
    for(int t = 0; t < 400; t++) {
        int n = 1 + rng() % 80, alphabet = 1 + rng() % 26;
        std::string s(n, 'a');
        for(char &c: s) c += rng() % alphabet;
        check(s);
    }
    check(std::string(300, 'a'));
    std::cout << "402 eertree cases passed brute-force palindrome, parent and suffix-link checks\n";
}
