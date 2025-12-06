#ifndef CP_ALGO_GEOMETRY_CONVEX_HULL_HPP
#define CP_ALGO_GEOMETRY_CONVEX_HULL_HPP
#include "point.hpp"
#include <algorithm>
#include <utility>
#include <vector>
#include <ranges>
namespace cp_algo::geometry{template<typename ftype>std::vector<point_t<ftype>>convex_hull(std::vector<point_t<ftype>>r){using point=point_t<ftype>;std::ranges::sort(r);if(size(r)<=1||r[0]==r.back()){return empty(r)?r:std::vector{r[0]};}std::vector<point>hull={r[0]};for(int half:{0,1}){size_t base=size(hull);for(auto it:std::views::drop(r,1)){while(size(hull)>=base+1){point a=hull.back();if(point::ccw(it-a,end(hull)[-2]-a)){break;}else{hull.pop_back();}}hull.push_back(it);}std::ranges::reverse(r);std::ignore=half;}hull.pop_back();return hull;}}
#endif
