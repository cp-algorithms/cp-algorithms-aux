#ifndef CP_ALGO_NUMBER_THEORY_PRIMALITY_HPP
#define CP_ALGO_NUMBER_THEORY_PRIMALITY_HPP
#include "modint.hpp"
#include <algorithm>
#include <bit>
namespace cp_algo::math{template<typename _Int>
bool is_prime(_Int m){using Int=std::make_signed_t<_Int>;
using UInt=std::make_unsigned_t<Int>;
if(m==1||m%2==0){return m==2;}
int s=std::countr_zero(UInt(m-1));
auto d=(m-1)>>s;
using base=dynamic_modint<Int>;
auto test=[&](base x){x=bpow(x,d);
if(std::abs(x.rem())<=1){return true;}
for(int i=1;i<s&&x!=-1;i++){x*=x;}
return x==-1;};
return base::with_mod(m,[&](){#ifdef CP_ALGO_NUMBER_THEORY_PRIMALITY_BASES_HPP
uint16_t base2=7,base3=61;
if(m!=uint32_t(m)){base2=base_table1[uint32_t(m*0xAD625B89)>>18];
base3=base_table2[base2>>13];}
return test(2)&&test(base2)&&test(base3);
#else
return std::ranges::all_of(std::array{2,325,9375,28178,450775,9780504,1795265022},test);});}}
#endif
#endif