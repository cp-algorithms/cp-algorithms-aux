#include <cassert>
#include <bits/stdc++.h>
#ifdef TEST_MINIFIED
#include "cp-algo/min/structures/treap.hpp"
#else
#include "cp-algo/structures/treap.hpp"
#endif
using node = cp_algo::structures::treap::node<cp_algo::structures::treap::null_meta>;
int main() {
    std::mt19937 rng(2741);
    for(int test = 0; test < 1000; ++test) {
        std::vector<node *> nodes;
        for(int i = 0, n = test == 0 ? 2 : 1 + rng() % 200; i < n; ++i) {
            nodes.push_back(new node({}, test == 0 ? 7 : rng() % 8));
        }
        auto root = node::build(nodes);
        assert(root->size == nodes.size());
        std::vector<node *> actual;
        node::exec_on_each(root, [&](node *p) { actual.push_back(p); });
        assert(actual == nodes);
        for(auto p : nodes) {
            for(auto child : p->children) {
                assert(!child || p->prior <= child->prior);
            }
        }
        auto [left, right] = node::split(root, nodes.size() / 2);
        root = node::merge(left, right);
        actual.clear();
        node::exec_on_each(root, [&](node *p) { actual.push_back(p); });
        assert(actual == nodes);
        for(auto p : nodes)
            delete p;
    }
}
