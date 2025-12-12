#ifndef CP_ALGO_STRUCTURES_EERTREE_HPP
#define CP_ALGO_STRUCTURES_EERTREE_HPP
#include "../util/big_alloc.hpp"
#include <forward_list>
#include <functional>
#include <iostream>
#include <vector>
#include <string>
namespace cp_algo::structures {
    template<int sigma = 26, char mch = 'a'>
    struct eertree {
        eertree(size_t q) {
            q += 2;
            s = std::string(q, -1);
            len = par = link = big_vector(q, 0);
            to.resize(q);
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
            if(!get(last, c)) {
                int u = get(get_link(link[last]), c);
                link[sz] = u;
                par[sz] = last;
                len[sz] = len[last] + 2;
                to[last].emplace_front((sz++ << 8) | c);
            }
            last = get(last, c);
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
        big_vector<std::forward_list<int>> to;
        big_vector<int> len, link, par;
        std::string s;
        int n = 1, sz = 2, last = 0;
    };
}
#endif // CP_ALGO_STRUCTURES_EERTREE_HPP
