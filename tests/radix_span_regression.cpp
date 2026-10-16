#include <cassert>
#include <bits/stdc++.h>
#ifdef TEST_MINIFIED
#include "cp-algo/min/util/compress_coords.hpp"
#else
#include "cp-algo/util/compress_coords.hpp"
#endif
struct item {
    int key, index;
    bool operator==(item const &) const = default;
};
int main() {
    std::vector<int> values = {1, 0};
    std::span<int> view(values);
    cp_algo::radix_sort(view);
    assert((values == std::vector<int>{0, 1}));
    values = {1, 0};
    std::vector<std::reference_wrapper<int>> refs(values.begin(), values.end());
    std::span<std::reference_wrapper<int>> references(refs);
    auto original = cp_algo::compress_coords(references);
    assert((std::vector<int>(original.begin(), original.end()) == std::vector<int>{0, 1}));
    assert((values == std::vector<int>{1, 0}));
    std::array<int, 4> counts = {3, 0, 2, 1};
    std::span<int> counted(counts);
    cp_algo::count_sort<4>(counted);
    assert((counts == std::array<int, 4>{0, 1, 2, 3}));
    std::mt19937 rng(3182);
    for(int test = 0; test < 500; ++test) {
        std::vector<item> items;
        for(int i = 0, n = 1 + rng() % 200; i < n; ++i)
            items.push_back({int(rng() % 9) - 4, i});
        auto expected = items;
        std::stable_sort(expected.begin(), expected.end(), [](auto a, auto b) { return a.key < b.key; });
        std::span<item> range(items);
        cp_algo::radix_sort(range, [](auto x) { return x.key; });
        assert(items == expected);
    }
}
