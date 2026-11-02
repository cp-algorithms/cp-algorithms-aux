#include <cassert>
#include <bits/stdc++.h>
#ifdef TEST_MINIFIED
#include "cp-algo/min/structures/fenwick_set.hpp"
#else
#include "cp-algo/structures/fenwick_set.hpp"
#endif
using set_type = cp_algo::structures::fenwick_set<128>;
void check(set_type const &values) {
    assert(values.sz == 2);
    assert(values.find_by_order(0) == 1);
    assert(values.find_by_order(1) == 70);
}
int main() {
    set_type original;
    original.insert(1);
    original.insert(70);
    auto copied = original;
    check(copied);
    auto const_copy = std::as_const(original);
    check(const_copy);
    auto moved = std::move(copied);
    check(moved);
    [](set_type value) { check(value); }(original);
    moved.erase(70);
    check(original);
    std::vector<int> keys = {1, 70};
    set_type from_range(keys);
    check(from_range);
}
