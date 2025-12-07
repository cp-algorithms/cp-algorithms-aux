#ifndef CP_ALGO_GEOMETRY_CONVEX_HULL_HPP
#define CP_ALGO_GEOMETRY_CONVEX_HULL_HPP
#include "point.hpp"
#include "../util/big_alloc.hpp"
#include <algorithm>
#include <utility>
#include <vector>
#include <ranges>
namespace cp_algo::geometry {
    auto convex_hull(auto r) {
        using point = std::decay_t<decltype(r[0])>;
        std::ranges::sort(r);
        if(size(r) <= 1 || r[0] == r.back()) {
            return empty(r) ? big_vector<point>{} : big_vector{r[0]};
        }
        big_vector<point> hull = {r[0]};
        for(int half: {0, 1}) {
            size_t base = size(hull);
            for(auto it: std::views::drop(r, 1)) {
                while(size(hull) >= base + 1) {
                    point a = hull.back();
                    if(point::ccw(it - a, end(hull)[-2] - a)) {
                        break;
                    } else {
                        hull.pop_back();
                    }
                }
                hull.push_back(it);
            }
            std::ranges::reverse(r);
            std::ignore = half;
        }
        hull.pop_back();
        return hull;
    }
}
#endif // CP_ALGO_GEOMETRY_CONVEX_HULL_HPP
