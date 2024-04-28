#ifndef CP_ALGO_GEOMETRY_POINT_HPP
#define CP_ALGO_GEOMETRY_POINT_HPP
#include <algorithm>
#include <iostream>
#include <complex>
#include <utility>
#include <vector>
#include <ranges>
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
        size_t n = size(r);
        if(n <= 1) {
            return r;
        }
        std::ranges::nth_element(r, std::begin(r));
        std::ranges::sort(std::views::drop(r, 1), point::ccw_abs, [&](auto it){return it - r[0];});
        if(r[0] == r.back()) {
            return {r[0]};
        }
        std::vector<point> res = {r[0]};
        for(size_t i = 1; i < n; i++) {
            while(res.size() > 1 && !point::ccw(res.back() - *(end(res) - 2), r[i] - res.back())) {
                res.pop_back();
            }
            res.push_back(r[i]);
        }
        return res;
    }
}
#endif // CP_ALGO_GEOMETRY_POINT_HPP
