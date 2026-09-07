#ifndef CP_ALGO_MATH_POLY_SPARSE_POW_HPP
#define CP_ALGO_MATH_POLY_SPARSE_POW_HPP
#include "impl.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math{template<typename T>poly_t<T>pow_sparse(poly_t<T>const&p,int64_t k,size_t n){assert(k>=0);if(n==0){return{};}if(k==0){return T(1);}if(p.is_zero()){return{};}size_t shift=p.trailing_xk();if(shift&&uint64_t(k)>(n-1)/shift){return{};}size_t offset=shift*uint64_t(k);auto q=poly::impl::pow_sparse_unit(p,T(k),n-offset,shift);q*=bpow(p.a[shift],k);q.mul_xk_inplace(offset);return q;}}
#pragma GCC pop_options
#endif
