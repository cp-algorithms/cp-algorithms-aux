#ifndef CP_ALGO_MATH_POLY_SERIES_LOG_HPP
#define CP_ALGO_MATH_POLY_SERIES_LOG_HPP
#include "inv.hpp"
#include "../calculus.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math{template<typename T>poly_t<T>log(poly_t<T>p,size_t n){if(n==0){return{};}assert(p[0]==T(1));p.mod_xk_inplace(n);auto dp=deriv(p);size_t k=n-1;if(k<magic){poly::impl::inv_inplace(p,k);p.mul_truncate(dp,k);return integr(std::move(p));}size_t m=std::bit_floor(k-1),t=k-m;auto r=inv(p,m);auto R=fft::dft<T>(r.a,m);typename poly_t<T>::Vector work(2*m);{auto D=fft::dft<T>(dp.a|std::views::take(m),m);D.mul(R,work,m);}poly_t<T>q(typename poly_t<T>::Vector(begin(work),begin(work)+m));{auto Q=fft::dft<T>(q.a,m);auto P=fft::dft<T>(p.a|std::views::take(k),m);P.mul_inplace(Q,work,k);}for(size_t i=0;i<t;i++){work[m+i]=dp[int(m+i)]-work[m+i];}auto E=fft::dft<T>(work|std::views::drop(m)|std::views::take(t),m);E.mul_inplace(R,work,t);q.a.resize(k);std::copy_n(begin(work),t,begin(q.a)+m);q.normalize();return integr(std::move(q));}}
#pragma GCC pop_options
#endif
