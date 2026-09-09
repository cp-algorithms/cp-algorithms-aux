#ifndef CP_ALGO_MATH_POLY_SERIES_EXP_HPP
#define CP_ALGO_MATH_POLY_SERIES_EXP_HPP
#include "inv.hpp"
#include "../calculus.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math{template<typename T>poly_t<T>exp(poly_t<T>p,size_t n){if(n==0){return{};}assert(p[0]==T(0));p.mod_xk_inplace(n);if(p.is_zero()){return T(1);}size_t m=std::min(n,std::bit_floor(size_t(magic-1)));typename poly_t<T>::Vector seed(m);seed[0]=1;for(size_t i=1;i<m;i++){for(size_t j=1;j<=i&&j<p.a.size();j++){seed[i]+=T(j)*p.a[j]*seed[i-j];}seed[i]*=small_inv<T>(i);}poly_t<T>q(std::move(seed));if(m==n){return q;}auto r=inv(q,m),dp=deriv(p);for(;m<n;m*=2){size_t k=std::min(2*m,n),t=k-m;auto Q=fft::dft<T>(q.a,m),R=fft::dft<T>(r.a,m);typename poly_t<T>::Vector work(2*m);{auto D=fft::dft<T>(dp.a|std::views::take(k-1),m);D.mul(Q,work,k-1);}auto E=fft::dft<T>(work|std::views::drop(m-1)|std::views::take(t),m);E.mul(R,work,t);for(size_t i=0;i<t;i++){work[i]*=small_inv<T>(m+i);}auto d=typename poly_t<T>::Vector(begin(work),begin(work)+t);{auto D=fft::dft<T>(d,m);D.mul(Q,work,t);}q.a.resize(k);std::copy_n(begin(work),t,begin(q.a)+m);if(k==n){break;}Q.mul(R,work,k);for(size_t i=0;i<t;i++){work[m+i]+=d[i];}auto H=fft::dft<T>(work|std::views::drop(m)|std::views::take(t),m);H.mul(R,work,t);r.a.resize(k);for(size_t i=0;i<t;i++){r.a[m+i]=-work[i];}}q.normalize();return q;}}
#pragma GCC pop_options
#endif
