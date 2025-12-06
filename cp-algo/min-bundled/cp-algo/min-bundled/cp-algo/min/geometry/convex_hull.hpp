#line 1 "cp-algo/min-bundled/cp-algo/min/geometry/convex_hull.hpp"
#line 1 "cp-algo/min/geometry/convex_hull.hpp"
#line 1 "cp-algo/min/geometry/point.hpp"
#line 1 "cp-algo/min/util/complex.hpp"
#include <iostream>
#include <cmath>
namespace cp_algo{template<typename T>
struct complex{using value_type=T;
T x,y;
constexpr complex():x(),y(){}
constexpr complex(T x):x(x),y(){}
constexpr complex(T x,T y):x(x),y(y){}
complex&operator*=(T t){x*=t;y*=t;return*this;}
complex&operator/=(T t){x/=t;y/=t;return*this;}
complex operator*(T t)const{return complex(*this)*=t;}
complex operator/(T t)const{return complex(*this)/=t;}
complex&operator+=(complex t){x+=t.x;y+=t.y;return*this;}
complex&operator-=(complex t){x-=t.x;y-=t.y;return*this;}
complex operator*(complex t)const{return{x*t.x-y*t.y,x*t.y+y*t.x};}
complex operator/(complex t)const{return*this*t.conj()/t.norm();}
complex operator+(complex t)const{return complex(*this)+=t;}
complex operator-(complex t)const{return complex(*this)-=t;}
complex&operator*=(complex t){return*this=*this*t;}
complex&operator/=(complex t){return*this=*this/t;}
complex operator-()const{return{-x,-y};}
complex conj()const{return{x,-y};}
T norm()const{return x*x+y*y;}
T abs()const{return std::sqrt(norm());}
T const real()const{return x;}
T const imag()const{return y;}
T&real(){return x;}
T&imag(){return y;}
static constexpr complex polar(T r,T theta){return{T(r*cos(theta)),T(r*sin(theta))};}
auto operator<=>(complex const&t)const=default;};
template<typename T>
complex<T>operator*(auto x,complex<T>y){return y*=x;}
template<typename T>complex<T>conj(complex<T>x){return x.conj();}
template<typename T>T norm(complex<T>x){return x.norm();}
template<typename T>T abs(complex<T>x){return x.abs();}
template<typename T>T&real(complex<T>&x){return x.real();}
template<typename T>T&imag(complex<T>&x){return x.imag();}
template<typename T>T const real(complex<T>const&x){return x.real();}
template<typename T>T const imag(complex<T>const&x){return x.imag();}
template<typename T>
constexpr complex<T>polar(T r,T theta){return complex<T>::polar(r,theta);}
template<typename T>
std::ostream&operator<<(std::ostream&out,complex<T>x){return out<<x.real()<<' '<<x.imag();}}
#line 5 "cp-algo/min/geometry/point.hpp"
namespace cp_algo::geometry{template<typename ftype>
struct point_t:complex<ftype>{using Base=complex<ftype>;
using Base::Base;
point_t(Base const&t):Base(t){}
auto operator<=>(point_t const&t)const{return std::pair{y(),-x()}<=>std::pair{t.y(),-t.x()};}
ftype x()const{return Base::real();}
ftype y()const{return Base::imag();}
point_t cmul(point_t const&t)const{return conj(*this)*t;}
ftype dot(point_t const&t)const{return cmul(t).x();}
ftype cross(point_t const&t)const{return cmul(t).y();}
static constexpr point_t O={0,0};
int half()const{return*this<O?-1:*this==O?0:1;}
static bool ccw(point_t const&a,point_t const&b){return a.cross(b)>0;}
static bool ccw_abs(point_t const&a,point_t const&b){return std::tuple{a.half(),(ftype)0,norm(a)}<
std::tuple{b.half(),a.cross(b),norm(b)};}
void read(){ftype _x,_y;
std::cin>>_x>>_y;
*this={_x,_y};}
void print()const{std::cout<<x()<<' '<<y()<<"\n";}};}
#line 4 "cp-algo/min/geometry/convex_hull.hpp"
#include <algorithm>
#include <utility>
#include <vector>
#include <ranges>
namespace cp_algo::geometry{template<typename ftype>
std::vector<point_t<ftype>>convex_hull(std::vector<point_t<ftype>>r){using point=point_t<ftype>;
std::ranges::sort(r);
if(size(r)<=1||r[0]==r.back()){return empty(r)?r:std::vector{r[0]};}
std::vector<point>hull={r[0]};
for(int half:{0,1}){size_t base=size(hull);
for(auto it:std::views::drop(r,1)){while(size(hull)>=base+1){point a=hull.back();
if(point::ccw(it-a,end(hull)[-2]-a)){break;}else{hull.pop_back();}}
hull.push_back(it);}
std::ranges::reverse(r);
std::ignore=half;}
hull.pop_back();
return hull;}}