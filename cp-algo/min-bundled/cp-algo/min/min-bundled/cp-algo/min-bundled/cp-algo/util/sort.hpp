#line 1 "cp-algo/min/min-bundled/cp-algo/min-bundled/cp-algo/util/sort.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/util/sort.hpp"
#line 1 "cp-algo/util/sort.hpp"
#line 1 "cp-algo/util/bit.hpp"
#line 1 "cp-algo/util/simd.hpp"
#include <experimental/simd>
#include <cstdint>
#include <cstddef>
#include <memory>
namespace cp_algo{template<typename T,size_t len>
using simd[[gnu::vector_size(len*sizeof(T))]]=T;
using i64x4=simd<int64_t,4>;
using u64x4=simd<uint64_t,4>;
using u32x8=simd<uint32_t,8>;
using i32x4=simd<int32_t,4>;
using u32x4=simd<uint32_t,4>;
using i16x4=simd<int16_t,4>;
using u8x32=simd<uint8_t,32>;
using dx4=simd<double,4>;
[[gnu::target("avx2")]]inline dx4 abs(dx4 a){return a<0?-a:a;}
static constexpr dx4 magic=dx4()+(3ULL<<51);
[[gnu::target("avx2")]]inline i64x4 lround(dx4 x){return i64x4(x+magic)-i64x4(magic);}
[[gnu::target("avx2")]]inline dx4 to_double(i64x4 x){return dx4(x+i64x4(magic))-magic;}
[[gnu::target("avx2")]]inline dx4 round(dx4 a){return dx4{std::nearbyint(a[0]),
std::nearbyint(a[1]),
std::nearbyint(a[2]),
std::nearbyint(a[3])};}
[[gnu::target("avx2")]]inline u64x4 low32(u64x4 x){return x&uint32_t(-1);}
[[gnu::target("avx2")]]inline auto swap_bytes(auto x){return decltype(x)(__builtin_shufflevector(u32x8(x),u32x8(x),1,0,3,2,5,4,7,6));}
[[gnu::target("avx2")]]inline u64x4 montgomery_reduce(u64x4 x,uint32_t mod,uint32_t imod){auto x_ninv=u64x4(_mm256_mul_epu32(__m256i(x),__m256i()+imod));
x+=u64x4(_mm256_mul_epu32(__m256i(x_ninv),__m256i()+mod));
return swap_bytes(x);}
[[gnu::target("avx2")]]inline u64x4 montgomery_mul(u64x4 x,u64x4 y,uint32_t mod,uint32_t imod){return montgomery_reduce(u64x4(_mm256_mul_epu32(__m256i(x),__m256i(y))),mod,imod);}
[[gnu::target("avx2")]]inline u32x8 montgomery_mul(u32x8 x,u32x8 y,uint32_t mod,uint32_t imod){return u32x8(montgomery_mul(u64x4(x),u64x4(y),mod,imod))|
u32x8(swap_bytes(montgomery_mul(u64x4(swap_bytes(x)),u64x4(swap_bytes(y)),mod,imod)));}
[[gnu::target("avx2")]]inline dx4 rotate_right(dx4 x){static constexpr u64x4 shuffler={3,0,1,2};
return __builtin_shuffle(x,shuffler);}
template<std::size_t Align=32>
[[gnu::target("avx2")]]inline bool is_aligned(const auto*p)noexcept{return(reinterpret_cast<std::uintptr_t>(p)%Align)==0;}
template<class Target>
[[gnu::target("avx2")]]inline Target&vector_cast(auto&&p){return*reinterpret_cast<Target*>(std::assume_aligned<alignof(Target)>(&p));}}
#line 5 "cp-algo/util/bit.hpp"
#include <array>
#include <bit>
namespace cp_algo{template<typename Uint>
constexpr size_t bit_width=sizeof(Uint)*8;
uint64_t mask(size_t n){return(1ULL<<n)-1;}
size_t order_of_bit(auto x,size_t k){return k?std::popcount(x<<(bit_width<decltype(x)>-k)):0;}
[[gnu::target("bmi2")]]inline size_t kth_set_bit(uint64_t x,size_t k){return std::countr_zero(_pdep_u64(1ULL<<k,x));}
template<int fl=0>
void with_bit_floor(size_t n,auto&&callback){if constexpr(fl>=63){return;}else if(n>>(fl+1)){with_bit_floor<fl+1>(n,callback);}else{callback.template operator()<1ULL<<fl>();}}
void with_bit_ceil(size_t n,auto&&callback){with_bit_floor(n,[&]<size_t N>(){if(N==n){callback.template operator()<N>();}else{callback.template operator()<N<<1>();}});}
[[gnu::target("avx2")]]inline uint32_t read_bits(char const*p){return _mm256_movemask_epi8(__m256i(vector_cast<u8x32 const>(p[0])+(127-'0')));}
[[gnu::target("avx2")]]inline uint64_t read_bits64(char const*p){return read_bits(p)|(uint64_t(read_bits(p+32))<<32);}
[[gnu::target("avx2")]]inline void write_bits(char*p,uint32_t bits){static constexpr u8x32 shuffler={0,0,0,0,0,0,0,0,
1,1,1,1,1,1,1,1,
2,2,2,2,2,2,2,2,
3,3,3,3,3,3,3,3};
auto shuffled=u8x32(_mm256_shuffle_epi8(__m256i()+bits,__m256i(shuffler)));
static constexpr u8x32 mask={1,2,4,8,16,32,64,128,
1,2,4,8,16,32,64,128,
1,2,4,8,16,32,64,128,
1,2,4,8,16,32,64,128};
for(int z=0;z<32;z++){p[z]=shuffled[z]&mask[z]?'1':'0';}}
[[gnu::target("avx2")]]inline void write_bits64(char*p,uint64_t bits){write_bits(p,uint32_t(bits));
write_bits(p+32,uint32_t(bits>>32));}}
#line 4 "cp-algo/util/sort.hpp"
#include <algorithm>
#include <numeric>
#include <ranges>
#include <vector>
namespace cp_algo{template<size_t maxc>
void count_sort(auto&a,auto&&proj){std::array<int,maxc>cnt={};
for(auto&x:a){cnt[proj(x)]++;}
std::partial_sum(begin(cnt),end(cnt),begin(cnt));
auto res=a;
for(auto const&it:a|std::views::reverse){res[--cnt[proj(it)]]=it;}
a=std::move(res);}
template<size_t maxc>
void count_sort(auto&a){count_sort<maxc>(a,std::identity{});}
void radix_sort(auto&a,auto&&proj){if(empty(a)){return;}
auto[mn,mx]=std::ranges::minmax(a,{},proj);
with_bit_floor<1>(size(a),[&]<size_t floor>(){constexpr int base=std::min<size_t>(floor,1<<16);
for(int64_t i=1;i<=std::invoke(proj,mx)-std::invoke(proj,mn);i*=base){count_sort<base>(a,[&](auto const&x){return(std::invoke(proj,x)-std::invoke(proj,mn))/i%base;});}});}
void radix_sort(auto&a){radix_sort(a,std::identity{});}}