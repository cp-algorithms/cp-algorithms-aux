#line 1 "cp-algo/min/min-bundled/cp-algo/min-bundled/cp-algo/util/simd.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/util/simd.hpp"
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