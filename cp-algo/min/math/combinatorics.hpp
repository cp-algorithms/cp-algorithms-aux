#ifndef CP_ALGO_MATH_COMBINATORICS_HPP
#define CP_ALGO_MATH_COMBINATORICS_HPP
#include "../math/common.hpp"
#include <cassert>
#include <ranges>
namespace cp_algo::math{template<typename T>
T fact(auto n){static std::vector<T>F(maxn);
static bool init=false;
if(!init){F[0]=T(1);
for(int i=1;i<maxn;i++){F[i]=F[i-1]*T(i);}
init=true;}
return F[n];}
template<typename T>
T rfact(auto n){static std::vector<T>F(maxn);
static bool init=false;
if(!init){int t=(int)std::min<int64_t>(T::mod(),maxn)-1;
F[t]=T(1)/fact<T>(t);
for(int i=t-1;i>=0;i--){F[i]=F[i+1]*T(i+1);}
init=true;}
return F[n];}
template<typename T,int base>
T pow_fixed(int n){static std::vector<T>prec_low(1<<16);
static std::vector<T>prec_high(1<<16);
static bool init=false;
if(!init){init=true;
prec_low[0]=prec_high[0]=T(1);
T step_low=T(base);
T step_high=bpow(T(base),1<<16);
for(int i=1;i<(1<<16);i++){prec_low[i]=prec_low[i-1]*step_low;
prec_high[i]=prec_high[i-1]*step_high;}}
return prec_low[n&0xFFFF]*prec_high[n>>16];}
template<typename T>
std::vector<T>bulk_invs(auto const&args){std::vector<T>res(std::size(args),args[0]);
for(size_t i=1;i<std::size(args);i++){res[i]=res[i-1]*args[i];}
auto all_invs=T(1)/res.back();
for(size_t i=std::size(args)-1;i>0;i--){res[i]=all_invs*res[i-1];
all_invs*=args[i];}
res[0]=all_invs;
return res;}
template<typename T>
T small_inv(auto n){static auto F=bulk_invs<T>(std::views::iota(1,maxn));
return F[n-1];}
template<typename T>
T binom_large(T n,auto r){assert(r<maxn);
T ans=1;
for(decltype(r)i=0;i<r;i++){ans=ans*T(n-i)*small_inv<T>(i+1);}
return ans;}
template<typename T>
T binom(auto n,auto r){if(r<0||r>n){return T(0);}else if(n>=maxn){return binom_large(T(n),r);}else{return fact<T>(n)*rfact<T>(r)*rfact<T>(n-r);}}}
#endif