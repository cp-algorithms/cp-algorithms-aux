#include <bits/stdc++.h>
#ifdef TEST_MINIFIED
#include "cp-algo/min/geometry/closest_pair.hpp"
#else
#include "cp-algo/geometry/closest_pair.hpp"
#endif
using point = cp_algo::geometry::point_t<int64_t>;
void check(std::vector<point> const& points) {
    auto [a, b] = cp_algo::geometry::closest_pair(points);
    int64_t expected = std::numeric_limits<int64_t>::max();
    for(size_t i = 0; i < points.size(); ++i) {
        for(size_t j = 0; j < i; ++j) {
            expected = std::min(expected, norm(points[i] - points[j]));
        }
    }
    assert(a != b && norm(points[a] - points[b]) == expected);
}
int main() {
    check({{-1000000000, -1000000000}, {1000000000, 1000000000}, {-1000000000, 1000000000}});
    check({{3, 4}, {3, 4}});
    std::mt19937_64 rng(371);
    for(int test = 0; test < 1000; ++test) {
        std::vector<point> points;
        for(int i = 0, n = 2 + rng() % 150; i < n; ++i) {
            points.emplace_back(int64_t(rng() % 2000000001) - 1000000000,
                                int64_t(rng() % 2000000001) - 1000000000);
        }
        check(points);
    }
}
