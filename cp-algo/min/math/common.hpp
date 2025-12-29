#ifndef CP_ALGO_MATH_COMMON_HPP
#define CP_ALGO_MATH_COMMON_HPP
#include <functional>
#include <cstdint>
#include <cassert>
#include <bit>
namespace cp_algo::math{
#ifdef CP_ALGO_MAXN
const int maxn=CP_ALGO_MAXN;
#else
const int maxn=1<<19;
#endif
const int magic=64;auto bpow(auto const&x,auto n,auto const&one,auto op){if(n==0){return one;}auto ans=x;for(int j=std::bit_width<uint64_t>(n)-2;~j;j--){ans=op(ans,ans);if((n>>j)&1){ans=op(ans,x);}}return ans;}auto bpow(auto x,auto n,auto ans){return bpow(x,n,ans,std::multiplies{});}template<typename T>T bpow(T const&x,auto n){return bpow(x,n,T(1));}inline constexpr auto inv2(auto x){assert(x%2);std::make_unsigned_t<decltype(x)>y=1;while(y*x!=1){y*=2-x*y;}return y;}}
#endif
