#ifndef CP_ALGO_MATH_POLY_POWMOD_HPP
#define CP_ALGO_MATH_POLY_POWMOD_HPP
#include "div.hpp"
#include "series.hpp"
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math{template<typename T>poly_t<T>circular_closure(poly_t<T>p,size_t m){assert(m>0);for(size_t i=p.a.size();i>m;--i){p.a[i-1-m]+=p.a[i-1];}p.mod_xk_inplace(m);return p;}template<typename T>poly_t<T>powmod_circular(poly_t<T>p,int64_t k,size_t m){assert(k>=0&&m>0);p=circular_closure(std::move(p),m);return bpow(p,k,poly_t<T>(1),[m](auto const&a,auto const&b){auto product=a;product*=&a==&b?product:b;return circular_closure(std::move(product),m);});}template<typename T>poly_t<T>powmod(poly_t<T>p,int64_t k,poly_t<T>const&md){assert(k>=0&&!md.is_zero());int d=md.deg();if(d==0){return{};}if(md==poly_t<T>::xk(d)){return pow(std::move(p),k,d);}if(md==poly_t<T>::xk(d)-poly_t<T>(1)){return powmod_circular(std::move(p),k,d);}auto mdri=inv(md.reversed(),d+1);return bpow(p%md,k,poly_t<T>(1),[&](auto const&a,auto const&b){auto product=a;product*=&a==&b?product:b;auto[q,r]=poly::impl::divmod_hint(std::move(product),md,mdri);return r;});}}
#pragma GCC pop_options
#endif
