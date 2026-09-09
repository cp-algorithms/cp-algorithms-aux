#ifndef CP_ALGO_MATH_POLY_SPARSE_INV_HPP
#define CP_ALGO_MATH_POLY_SPARSE_INV_HPP
#include "impl.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math{template<typename T>poly_t<T>inv_sparse(poly_t<T>const&p,size_t n){if(n==0){return{};}assert(p[0]!=T(0));auto terms=poly::impl::sparse_terms(p,n);T a0inv=T(1)/p[0];for(auto&[j,a]:terms){a*=-a0inv;}typename poly_t<T>::Vector q(n);q[0]=a0inv;for(size_t i=1;i<n;i++){q[i]=poly::impl::sparse_dot<T>(terms,q,i);}return q;}}
#pragma GCC pop_options
#endif
