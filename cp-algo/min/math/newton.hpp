#ifndef CP_ALGO_MATH_NEWTON_HPP
#define CP_ALGO_MATH_NEWTON_HPP
#include "poly.hpp"
namespace cp_algo::math{template<typename base>auto newton(auto&&FFd,base f0,size_t n){using polyn=poly_t<base>;polyn f=polyn(f0);for(size_t len=1;len<n;len*=2){auto[Ff,Fdf]=FFd(f,2*len);Ff.div_xk_inplace(len);Ff.mul_truncate(Fdf.inv_inplace(len),len);f-=Ff.mul_xk_inplace(len);}return f.mod_xk_inplace(n);}}
#endif
