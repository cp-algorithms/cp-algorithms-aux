#ifndef CP_ALGO_MATH_POLY_SPARSE_SQRT_HPP
#define CP_ALGO_MATH_POLY_SPARSE_SQRT_HPP
#include "impl.hpp"
#include "../../../number_theory/discrete_sqrt.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math{template<typename T>std::optional<poly_t<T>>sqrt_sparse(poly_t<T>const&p,size_t n){if(n==0||p.is_zero()){return poly_t<T>{};}size_t shift=p.trailing_xk();if(shift>=n){return poly_t<T>{};}if(shift%2){return std::nullopt;}auto c=math::sqrt(p.a[shift]);if(!c){return std::nullopt;}auto q=poly::impl::pow_sparse_unit(p,T(1)/T(2),n-shift,shift);q*=*c;q.mul_xk_inplace(shift/2);return q;}}
#pragma GCC pop_options
#endif
