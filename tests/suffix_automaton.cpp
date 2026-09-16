#include "cp-algo/structures/suffix_automaton.hpp"
#include <iostream>
#include <numeric>
#include <random>
#include <set>
using cp_algo::structures::suffix_automaton;

void check_match(suffix_automaton const &sam, std::string const &s, std::string const &t, int best) {
    auto [a, b, c, d] = sam.longest_common_substring(t);
    assert(0 <= a && a <= b && b <= (int)s.size());
    assert(0 <= c && c <= d && d <= (int)t.size());
    assert(b - a == best && d - c == best);
    assert(s.substr(a, b - a) == t.substr(c, d - c));
}
void check_suffix_order(std::string const &s) {
    std::vector<int> expected(s.size());
    std::iota(expected.begin(), expected.end(), 0);
    std::string_view text(s);
    std::ranges::sort(expected, [&](int i, int j) { return text.substr(i) < text.substr(j); });
    assert(std::ranges::equal(cp_algo::structures::suffix_array(s), expected));
}
void check(std::string const &s, std::string const &t) {
    suffix_automaton sam(s);
    std::set<std::string> substrings;
    for(size_t i = 0; i < s.size(); i++) {
        for(size_t j = i + 1; j <= s.size(); j++) substrings.insert(s.substr(i, j - i));
    }
    assert(sam.count_distinct() == (int64_t)substrings.size());
    int best = 0;
    for(auto const &part: substrings) if(t.find(part) != std::string::npos) best = std::max(best, (int)part.size());
    check_match(sam, s, t, best);
    check_match(sam, s, s, (int)s.size());
    check_match(sam, s, "#", 0);
    assert(sam.count_distinct() == (int64_t)substrings.size());
    auto copy = sam;
    check_match(copy, s, t, best);
    check_suffix_order(s);
}
int64_t reference_count(std::string const &s) {
    struct node { int len = 0, link = -1; std::array<int, 26> next{}; };
    std::vector<node> states(1);
    states.reserve(2 * s.size() + 1);
    int last = 0;
    for(char c: s) {
        int x = c - 'a', current = (int)states.size();
        states.push_back({states[last].len + 1, 0, {}});
        int p = last;
        while(p >= 0 && !states[p].next[x]) { states[p].next[x] = current; p = states[p].link; }
        if(p >= 0) {
            int q = states[p].next[x];
            if(states[p].len + 1 == states[q].len) states[current].link = q;
            else {
                int clone = (int)states.size();
                auto row = states[q]; row.len = states[p].len + 1;
                states.push_back(row);
                while(p >= 0 && states[p].next[x] == q) { states[p].next[x] = clone; p = states[p].link; }
                states[q].link = states[current].link = clone;
            }
        }
        last = current;
    }
    int64_t answer = 0;
    for(size_t i = 1; i < states.size(); i++) answer += states[i].len - states[states[i].link].len;
    return answer;
}
int main() {
    std::mt19937 rng(831);
    for(int len = 0; len <= 10; len++) {
        for(int mask = 0; mask < (1 << len); mask++) {
            std::string s(len, 'a');
            for(int i = 0; i < len; i++) s[i] += (mask >> i) & 1;
            check(s, std::string(s.rbegin(), s.rend()));
        }
    }
    for(int alphabet: {1, 2, 3, 4, 5, 8, 16, 17, 26}) {
        for(int rep = 0; rep < 30; rep++) {
            std::string s(rng() % 65, 'a'), t(rng() % 80, 'a');
            for(char &c: s) c += rng() % alphabet;
            for(char &c: t) c += rng() % alphabet;
            check(s, t);
        }
    }
    // Long clones exercise the cached-length fallback; borders exercise skipped leaf ranges.
    for(int n: {254, 255, 256, 511, 1024}) {
        for(auto const &s: {std::string(n, 'a'), "a" + std::string(n, 'b'),
                            "ab" + std::string(n, 'a'), std::string(n, 'a') + "b" + std::string(n, 'a')})
            check_suffix_order(s);
        std::string periodic;
        for(int i = 0; i < n; i++) periodic += char('a' + i % 3);
        check_suffix_order(periodic);
    }
    for(int alphabet: {17, 26}) {
        int boundary = (1 << 23) / alphabet;
        for(int n: {boundary - 1, boundary + 1}) {
            std::string s(n, 'a');
            for(char &c: s) c += rng() % alphabet;
            suffix_automaton sam(s);
            auto expected = reference_count(s);
            assert(sam.count_distinct() == expected);
            check_match(sam, s, s, n);
            check_match(sam, s, s.substr(10, 200) + '#' + s.substr(n - 60), 200);
            assert(sam.count_distinct() == expected);
        }
    }
    std::cout << "suffix automaton tests passed\n";
}
