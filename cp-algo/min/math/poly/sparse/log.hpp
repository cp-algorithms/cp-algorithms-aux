#ifndef CP_ALGO_MATH_POLY_SPARSE_LOG_HPP
#define CP_ALGO_MATH_POLY_SPARSE_LOG_HPP
#include "impl.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math{template<typename T>poly_t<T>log_sparse(poly_t<T>const&p,size_t n){if(n==0){return{};}assert(p[0]==T(1));auto terms=poly::impl::sparse_terms(p,n);typename poly_t<T>::Vector q(n);for(size_t i=1;i<n;i++){q[i]=T(i)*p[int(i)];q[i]-=poly::impl::sparse_dot<T>(terms,q,i);}for(size_t i=1;i<n;i++){q[i]*=small_inv<T>(i);}return q;}}
#pragma GCC pop_options
#endif
