#ifndef CP_ALGO_STRUCTURES_BIT_ARRAY_UTIL_HPP
#define CP_ALGO_STRUCTURES_BIT_ARRAY_UTIL_HPP
#include "../structures/bit_array.hpp"
#include "../util/bit.hpp"
#include <cstdint>
#include <cstddef>
#include <string>
#include <array>
#include <span>
CP_ALGO_BIT_PRAGMA_PUSHnamespace cp_algo::structures{template<typename BitArray>constexpr void from_string(BitArray&arr,std::span<char>bits){arr.resize(std::size(bits));int64_t i=0;constexpr int width=BitArray::width;for(;i+width<=(int64_t)std::size(bits);i+=width){arr.word(i/width)=read_bits64(bits.data()+i);}for(;i<(int64_t)std::size(bits);i++){if(bits[i]&1){arr.set(i);}}}template<typename BitArray>constexpr big_string to_string(BitArray const&arr){big_string res(arr.words*BitArray::width,'0');for(size_t i=0;i<arr.words;i++){write_bits64(res.data()+i*BitArray::width,arr.word(i));}res.resize(arr.n);return res;}template<typename BitArray>constexpr size_t count(BitArray const&arr,size_t n){size_t res=0;for(size_t i=0;i<n/BitArray::width;i++){res+=std::popcount(arr.word(i));}if(n%BitArray::width){res+=std::popcount(arr.word(n/BitArray::width)&mask(n%BitArray::width));}return res;}template<typename BitArray>constexpr size_t count(BitArray const&arr){return count(arr,arr.n);}template<typename BitArray>constexpr size_t ctz(BitArray const&arr){size_t res=0;size_t i=0;while(i<arr.words&&arr.word(i)==0){res+=BitArray::width;i++;}if(i<arr.words){res+=std::countr_zero(arr.word(i));}return std::min(res,arr.n);}template<typename BitArray>constexpr size_t skip(BitArray const&arr,size_t pos=0,int k=0){if(pos>=arr.n)return arr.n;size_t word_idx=pos/BitArray::width;auto w=arr.word(word_idx)>>(pos%BitArray::width);auto popcnt=std::popcount(w);if(popcnt>k){return pos+cp_algo::kth_set_bit(w,k);}k-=popcnt;while(++word_idx<arr.words){w=arr.word(word_idx);auto popcnt=std::popcount(w);if(popcnt>k)[[unlikely]]{return word_idx*BitArray::width+cp_algo::kth_set_bit(w,k);}k-=popcnt;}return arr.n;}}
#pragma GCC pop_options
#endif
