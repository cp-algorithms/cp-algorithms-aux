#ifndef CP_ALGO_MATH_POLY_SERIES_SQRT_HPP
#define CP_ALGO_MATH_POLY_SERIES_SQRT_HPP
#include "inv.hpp"
#include "../../../number_theory/discrete_sqrt.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math{template<typename T>std::optional<poly_t<T>>sqrt(poly_t<T>p,size_t n){if(n==0){return poly_t<T>{};}p.mod_xk_inplace(n);if(p.is_zero()){return p;}size_t shift=p.trailing_xk();if(shift%2){return std::nullopt;}if(shift){p.div_xk_inplace(shift);auto ans=sqrt(std::move(p),n-shift);if(ans){ans->mul_xk_inplace(shift/2);}return ans;}auto c=math::sqrt(p[0]);if(!c){return std::nullopt;}size_t m=std::min(n,std::bit_floor(size_t(magic-1)));typename poly_t<T>::Vector seed(m);T half_inv=T(1)/(T(2)**c);seed[0]=*c;for(size_t i=1;i<m;i++){seed[i]=p[int(i)];for(size_t j=1;j<i;j++){seed[i]-=seed[j]*seed[i-j];}seed[i]*=half_inv;}poly_t<T>ans(std::move(seed));if(m==n){return ans;}auto r=inv(ans,m);T half=T(1)/T(2);for(;m<n;m*=2){size_t k=std::min(2*m,n),t=k-m;auto R=fft::dft<T>(r.a,m);typename poly_t<T>::Vector work(2*m);{auto A=fft::dft<T>(ans.a,m);A.mul(A,work,k);}for(size_t i=0;i<t;i++){work[m+i]=(p[int(m+i)]-work[m+i])*half;}{auto E=fft::dft<T>(work|std::views::drop(m)|std::views::take(t),m);E.mul(R,work,t);}ans.a.resize(k);std::copy_n(begin(work),t,begin(ans.a)+m);if(k==n){break;}auto A=fft::dft<T>(ans.a,m);A.mul(R,work,k);auto E=fft::dft<T>(work|std::views::drop(m)|std::views::take(t),m);E.mul(R,work,t);r.a.resize(k);for(size_t i=0;i<t;i++){r.a[m+i]=-work[i];}}ans.normalize();return ans;}}
#pragma GCC pop_options
#endif
