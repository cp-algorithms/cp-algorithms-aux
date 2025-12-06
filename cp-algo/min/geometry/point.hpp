#ifndef CP_ALGO_GEOMETRY_POINT_HPP
#define CP_ALGO_GEOMETRY_POINT_HPP
#include "../util/complex.hpp"
#include <iostream>
namespace cp_algo::geometry{
template<typename ftype>
struct point_t:complex<ftype>{
using Base=complex<ftype>;
using Base::Base;
point_t(Base const&t):Base(t){}
auto operator<=>(point_t const&t)const{
return std::pair{y(),-x()}<=>std::pair{t.y(),-t.x()};
}
ftype x()const{return Base::real();}
ftype y()const{return Base::imag();}
point_t cmul(point_t const&t)const{return conj(*this)*t;}
ftype dot(point_t const&t)const{return cmul(t).x();}
ftype cross(point_t const&t)const{return cmul(t).y();}
static constexpr point_t O={0,0};
int half()const{
return*this<O?-1:*this==O?0:1;
}
static bool ccw(point_t const&a,point_t const&b){
return a.cross(b)>0;
}
static bool ccw_abs(point_t const&a,point_t const&b){
return std::tuple{a.half(),(ftype)0,norm(a)}<
std::tuple{b.half(),a.cross(b),norm(b)};
}
void read(){
ftype _x,_y;
std::cin>>_x>>_y;
*this={_x,_y};
}
void print()const{
std::cout<<x()<<' '<<y()<<"\n";
}
};
}
#endif
