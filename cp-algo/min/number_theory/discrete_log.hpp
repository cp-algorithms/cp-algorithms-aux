#ifndef CP_ALGO_NUMBER_THEORY_DISCRETE_LOG_HPP
#define CP_ALGO_NUMBER_THEORY_DISCRETE_LOG_HPP
#include "euler.hpp"
#include <optional>
namespace cp_algo::math{template<typename _Int>std::optional<_Int>discrete_log(_Int b,_Int c,_Int m,_Int a=1){if(std::abs(a-c)%m==0){return 0;}if(std::gcd(a,m)!=std::gcd(int64_t(a)*b,int64_t(m))){auto res=discrete_log(b,c,m,_Int(int64_t(a)*b%m));return res?std::optional(*res+1):res;}using Int=std::make_signed_t<_Int>;using base=dynamic_modint<Int>;return base::with_mod(m,[&]()->std::optional<_Int>{int sqrtmod=std::max(1,(int)std::sqrt(m)/2);big_map<_Int,int>small;base cur=a;for(int i=0;i<sqrtmod;i++){small[cur.getr()]=i;cur*=b;}base step=bpow(base(b),sqrtmod);cur=1;for(ptrdiff_t k=0;k<m;k+=sqrtmod){auto it=small.find((base(c)*cur).getr());if(it!=end(small)){auto cand=base::with_mod(period(base(b)),[&](){return base(it->second-k).getr();});if(base(a)*bpow(base(b),cand)==base(c)){return cand;}else{return std::nullopt;}}cur*=step;}return std::nullopt;});}}
#endif
