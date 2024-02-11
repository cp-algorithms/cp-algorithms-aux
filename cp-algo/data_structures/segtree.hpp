#ifndef CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_HPP
#define CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_HPP
#include <vector>
namespace cp_algo::data_structures {
    template<typename meta>
    struct segtree_t {
        const int N;
        std::vector<meta> _meta;

        segtree_t(int n): N(n), _meta(4 * N) {}

        segtree_t(std::vector<meta> leafs): N(size(leafs)), _meta(4 * N) {
            build(leafs);
        }

        void pull(int v, int l, int r) {
            if(r - l > 1) {
                _meta[v].pull(_meta[2 * v], _meta[2 * v + 1], l, r);
            }
        }

        void push(int v, int l, int r) {
            if(r - l > 1) {
                _meta[v].push(&_meta[2 * v], &_meta[2 * v + 1], l, r);
            } else {
                _meta[v].push(nullptr, nullptr, l, r);
            }
        }

        void build(auto &a, int v, size_t l, size_t r) {
            if(r - l == 1) {
                if(l < size(a)) {
                    _meta[v] = a[l];
                }
            } else {
                size_t m = (l + r) / 2;
                build(a, 2 * v, l, m);
                build(a, 2 * v + 1, m, r);
                pull(v, l, r);
            }
        }

        void build(auto &a) {
            build(a, 1, 0, N);
        }

        void exec_on_segment(int a, int b, auto func, auto proceed, auto stop, int v, int l, int r) {
            push(v, l, r);
            if(r <= a || b <= l || stop(_meta[v])) {
                return;
            } else if(a <= l && r <= b && proceed(_meta[v])) {
                func(_meta[v]);
                push(v, l, r);
            } else {
                int m = (l + r) / 2;
                exec_on_segment(a, b, func, proceed, stop, 2 * v, l, m);
                exec_on_segment(a, b, func, proceed, stop, 2 * v + 1, m, r);
                pull(v, l, r);
            }
        }

        static constexpr auto default_true = [](auto const&){return true;};
        static constexpr auto default_false = [](auto const&){return false;};

        void exec_on_segment(int a, int b, auto func, auto proceed, auto stop) {
            exec_on_segment(a, b, func, proceed, stop, 1, 0, N);
        }

        void exec_on_segment(int a, int b, auto func) {
            exec_on_segment(a, b, func, default_true, default_false);
        }
    };
}
#endif // CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_HPP