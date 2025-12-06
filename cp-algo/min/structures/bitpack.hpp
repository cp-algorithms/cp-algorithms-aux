#ifndef CP_ALGO_STRUCTURES_BITPACK_HPP
#define CP_ALGO_STRUCTURES_BITPACK_HPP
#include "../structures/bit_array.hpp"
#include "../util/simd.hpp"
#include <cstdint>
#include <cstddef>
#include <string>
#include <array>
namespace cp_algo::structures{template<typename BitArray>
struct _bitpack:BitArray{using Base=BitArray;
using Base::Base,Base::width,Base::words,Base::data,Base::n,Base::word;
auto operator<=>(_bitpack const&t)const=default;
constexpr _bitpack(std::string&bits):_bitpack(std::size(bits)){bits+=std::string(-std::size(bits)%width,'0');
for(size_t i=0;i<words;i++){word(i)=read_bits64(bits.data()+i*width);}}
constexpr _bitpack&xor_hint(_bitpack const&t,size_t hint){for(size_t i=hint/width;i<std::size(data);i++){data[i]^=t.data[i];}
return*this;}
constexpr _bitpack&operator^=(_bitpack const&t){return xor_hint(t,0);}
constexpr _bitpack operator^(_bitpack const&t)const{return _bitpack(*this)^=t;}
constexpr std::string to_string()const{std::string res(words*width,'0');
for(size_t i=0;i<words;i++){write_bits64(res.data()+i*width,word(i));}
res.resize(n);
return res;}
constexpr size_t count(size_t n)const{size_t res=0;
for(size_t i=0;i<n/width;i++){res+=std::popcount(word(i));}
if(n%width){res+=std::popcount(word(n/width)&mask(n%width));}
return res;}
constexpr size_t count()const{return count(n);}
constexpr size_t ctz()const{size_t res=0;
size_t i=0;
while(i<words&&word(i)==0){res+=width;
i++;}
if(i<words){res+=std::countr_zero(word(i));}
return std::min(res,n);}};
template<size_t N>
using bitpack=_bitpack<bit_array<N>>;
using dynamic_bitpack=_bitpack<dynamic_bit_array>;}
#endif