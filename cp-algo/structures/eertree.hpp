#ifndef CP_ALGO_STRUCTURES_EERTREE_HPP
#define CP_ALGO_STRUCTURES_EERTREE_HPP
#include "../util/big_alloc.hpp"
#include "stack_union.hpp"
#include <array>
#include <functional>
#include <iostream>
#include <vector>
#include <string>
namespace cp_algo::structures {
    template<int sigma = 26, char mch = 'a'>
    struct eertree {
        eertree(size_t q) {
            q += 2;
            s = big_string(q, -1);
            len = par = link = big_vector(q, 0);
            to = stack_union<int>((int)q);
            to.reserve((int)q);
            link[0] = 1;
            len[1] = -1;
        }
        
        int get_link(int v) const {
            while(s[n - 1] != s[n - len[v] - 2]) {
                v = link[v];
            }
            return v;
        }
        
        int get(int v, int c) const {
            if(v < 2 && c == char(c)) return root_to[v][(unsigned char)c];
            for(int cu: to[v]) {
                if(char(cu) == c) {
                    return cu >> 8;
                }
            }
            return 0;
        }
        
        void add_letter(char c) {
            c -= 'a';
            s[n++] = c;
            last = get_link(last);
            int v = get(last, c);
            if(!v) {
                v = sz++;
                link[v] = get(get_link(link[last]), c);
                par[v] = last;
                len[v] = len[last] + 2;
                to.push(last, (v << 8) | c);
                if(last < 2) root_to[last][(unsigned char)c] = v;
            }
            last = v;
        }
        int sufpal(auto &&adjust) const {
            return adjust(last);
        }
        int sufpal() const {
            return sufpal(std::identity{});
        }
        void print(auto &&adjust) const {
            std::cout << sz - 2 << "\n";
            for(int i = 2; i < sz; i++) {
                std::cout << adjust(par[i]) << ' ' << adjust(link[i]) << "\n";
            }
        }
        void print() const {
            print(std::identity{});
        }
    private:
        // Cache the two hot roots; other transitions stay in compact lists.
        std::array<std::array<int, 256>, 2> root_to{};
        stack_union<int> to;
        big_vector<int> len, link, par;
        big_string s;
        int n = 1, sz = 2, last = 0;
    };
}
#endif // CP_ALGO_STRUCTURES_EERTREE_HPP
