#ifndef CP_ALGO_GEOMETRY_POINT_HPP
#define CP_ALGO_GEOMETRY_POINT_HPP
#include "../random/rng.hpp"
#include <algorithm>
#include <iostream>
#include <complex>
#include <utility>
#include <vector>
#include <ranges>
#include <map>
namespace cp_algo::geometry {
    template<typename ftype>
    struct point_t: public std::complex<ftype> {
        using Base = std::complex<ftype>;
        using Base::Base;

        point_t(Base const& t): Base(t) {}
        
        auto operator <=> (point_t const& t) const {
            return std::pair{y(), -x()} <=> std::pair{t.y(), -t.x()};
        }

        ftype x() const {return Base::real();}
        ftype y() const {return Base::imag();}

        point_t cmul(point_t const& t) const {return conj(*this) * t;}
        ftype dot(point_t const& t) const {return cmul(t).x();}
        ftype cross(point_t const& t) const {return cmul(t).y();}

        static constexpr point_t O = {0, 0};

        int half() const {
            return *this < O ? -1 : *this == O ? 0 : 1;
        }

        static bool ccw(point_t const& a, point_t const& b) {
            return a.cross(b) > 0;
        }
        static bool ccw_abs(point_t const& a, point_t const & b) {
            return std::tuple{a.half(), (ftype)0, norm(a)} <
                   std::tuple{b.half(), a.cross(b), norm(b)};
        }
        void read() {
            ftype _x, _y;
            std::cin >> _x >> _y;
            *this = {_x, _y};
        }
        void print() const {
            std::cout << x() << ' ' << y() << "\n";
        }
    };

    template<typename point>
    std::vector<point> convex_hull(std::vector<point> r) {
        std::ranges::sort(r);
        if(size(r) <= 1 || r[0] == r.back()) {
            return empty(r) ? r : std::vector{r[0]};
        }
        std::vector<point> hull = {r[0]};
        for(int half: {0, 1}) {
            size_t base = size(hull);
            for(auto it: std::views::drop(r, 1)) {
                while(size(hull) >= base + 1) {
                    point a = hull.back();
                    if(point::ccw(it - a, *(end(hull) - 2) - a)) {
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

    // Rabin & Lipton
    template<typename point>
    auto closest_pair(std::vector<point> const& r) {
        size_t n = size(r);
        int64_t md = 1e18;
        for(size_t i = 0; i < n; i++) {
            auto A = random::rng() % n;
            auto B = random::rng() % n;
            if(A != B) {
                md = std::min(md, norm(r[A] - r[B]));
                if(md == 0) {
                    return std::pair{A, B};
                }
            }
        }
        std::map<point, std::vector<int>> neigs;
        md = ceil(sqrtl(md));
        for(size_t i = 0; i < n; i++) {
            neigs[r[i] / md].push_back(i);
        }
        size_t a = 0, b = 1;
        md = norm(r[a] - r[b]);
        for(auto &[p, id]: neigs) {
            for(int dx: {-1, 0, 1}) {
                for(int dy: {-1, 0, 1}) {
                    auto pp = p + point{dx, dy};
                    if(!neigs.count(pp)) {
                        continue;
                    }
                    for(size_t i: neigs[p + point{dx, dy}]) {
                        for(size_t j: id) {
                            if(j == i) {
                                break;
                            }
                            int64_t cur = norm(r[i] - r[j]);
                            if(cur < md) {
                                md = cur;
                                a = i;
                                b = j;
                            }
                        }
                    }
                }
            }
        }
        return std::pair{a, b};
    }
}
#endif // CP_ALGO_GEOMETRY_POINT_HPP
