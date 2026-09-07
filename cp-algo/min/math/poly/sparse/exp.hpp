#ifndef CP_ALGO_MATH_POLY_SPARSE_EXP_HPP
#define CP_ALGO_MATH_POLY_SPARSE_EXP_HPP
#include "impl.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math{template<typename T>poly_t<T>exp_sparse(poly_t<T>const&p,size_t n){if(n==0){return{};}assert(p[0]==T(0));auto terms=poly::impl::sparse_terms(p,n);for(auto&[j,a]:terms){a*=T(j);}typename poly_t<T>::Vector q(n);q[0]=1;for(size_t i=1;i<n;i++){q[i]=poly::impl::sparse_dot<T>(terms,q,i);q[i]*=small_inv<T>(i);}return q;}}
#pragma GCC pop_options
#endif
