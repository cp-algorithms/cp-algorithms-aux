#ifndef CP_ALGO_MATH_POLY_SERIES_LOG_HPP
#define CP_ALGO_MATH_POLY_SERIES_LOG_HPP
#include "inv.hpp"
#include "../calculus.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math{template<typename T>poly_t<T>log(poly_t<T>p,size_t n){if(n==0){return{};}assert(p[0]==T(1));p.mod_xk_inplace(n);auto dp=deriv(p);size_t k=n-1;if(k<magic){poly::impl::inv_inplace(p,k);p.mul_truncate(dp,k);return integr(std::move(p));}size_t m=std::bit_floor(k-1),t=k-m;auto r=inv(p,m);auto R=fft::spectrum<T>(r.a,2*m);typename poly_t<T>::Vector work(2*m);fft::spectrum<T>(dp.a|std::views::take(m),2*m).multiply(R,work,m);poly_t<T>q(typename poly_t<T>::Vector(begin(work),begin(work)+m));{auto Q=fft::spectrum<T>(q.a,2*m);fft::spectrum<T>(p.a|std::views::take(k),2*m).multiply(Q,work,k);}for(size_t i=0;i<t;i++){work[m+i]=dp[int(m+i)]-work[m+i];}fft::spectrum<T>(work|std::views::drop(m)|std::views::take(t),2*m).multiply(R,work,t);q.a.resize(k);std::copy_n(begin(work),t,begin(q.a)+m);q.normalize();return integr(std::move(q));}}
#pragma GCC pop_options
#endif
