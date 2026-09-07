#ifndef CP_ALGO_MATH_POLY_DIV_HPP
#define CP_ALGO_MATH_POLY_DIV_HPP
#include "impl/div.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math{template<typename T>std::array<poly_t<T>,2>divmod(poly_t<T>p,poly_t<T>const&q){assert(!q.is_zero());int n=p.deg()-q.deg();if(std::min(n,q.deg())<magic){return poly::impl::divmod_slow(std::move(p),q);}auto qi=inv(q.reversed(),n+1);return poly::impl::divmod_hint(std::move(p),q,qi);}}
#pragma GCC pop_options
#endif
