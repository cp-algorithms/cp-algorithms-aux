#line 1 "verify/structures/bitpack/prod_mod_2.test.cpp"
#define PROBLEM "https:#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_CHECKPOINT
#line 1 "cp-algo/structures/bitpack.hpp"
#line 1 "cp-algo/structures/bit_array.hpp"
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
#line 4 "cp-algo/structures/bit_array.hpp"
#include <cassert>
namespace cp_algo::structures{template<typename C>
concept Resizable=requires(C&c,std::size_t n){c.resize(n);};
template<class Cont>
struct _bit_array{static constexpr size_t width=bit_width<uint64_t>;
size_t words,n;
alignas(32)Cont data;
constexpr void resize(size_t N){n=N;
words=(n+width-1)/width;
if constexpr(Resizable<Cont>){data.resize(words);}else{assert(std::size(data)>=words);}}
constexpr _bit_array():n(0),words(0),data(){}
constexpr _bit_array(size_t N):data(){resize(N);}
constexpr uint64_t&word(size_t x){return data[x];}
constexpr uint64_t word(size_t x)const{return data[x];}
constexpr void set_all(uint64_t val=-1){for(auto&w:data){w=val;}}
constexpr void reset(){set_all(0);}
constexpr void set(size_t x){word(x/width)|=1ULL<<(x%width);}
constexpr void reset(size_t x){word(x/width)&=~(1ULL<<(x%width));}
constexpr void flip(size_t x){word(x/width)^=1ULL<<(x%width);}
constexpr bool test(size_t x)const{return(word(x/width)>>(x%width))&1;}
constexpr bool operator[](size_t x)const{return test(x);}
constexpr size_t size()const{return n;}};
template<size_t N>
struct bit_array:_bit_array<std::array<uint64_t,(N+63)/64>>{using Base=_bit_array<std::array<uint64_t,(N+63)/64>>;
using Base::Base,Base::words,Base::data;
constexpr bit_array():Base(N){}};
struct dynamic_bit_array:_bit_array<std::vector<uint64_t>>{using Base=_bit_array<std::vector<uint64_t>>;
using Base::Base,Base::words;
constexpr dynamic_bit_array(size_t N):Base(N){data.resize(words);}};}
#line 7 "cp-algo/structures/bitpack.hpp"
#include <string>
#line 9 "cp-algo/structures/bitpack.hpp"
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
#line 1 "cp-algo/util/checkpoint.hpp"
#include <iostream>
#include <chrono>
#line 6 "cp-algo/util/checkpoint.hpp"
#include <map>
namespace cp_algo{std::map<std::string,double>checkpoints;
template<bool final=false>
void checkpoint([[maybe_unused]]std::string const&msg=""){#ifdef CP_ALGO_CHECKPOINT
static double last=0;
double now=(double)clock()/CLOCKS_PER_SEC;
double delta=now-last;
last=now;
if(msg.size()&&!final){checkpoints[msg]+=delta;}
if(final){for(auto const&[key,value]:checkpoints){std::cerr<<key<<": "<<value*1000<<" ms\n";}
std::cerr<<"Total: "<<now*1000<<" ms\n";}
#endif}}
#line 7 "verify/structures/bitpack/prod_mod_2.test.cpp"
#include <bits/stdc++.h>
using namespace std;
const int maxn=1<<12;
const size_t K=8;
using bitpack=cp_algo::structures::bitpack<maxn>;
bitpack a[maxn],b[maxn],c[maxn];
bitpack precalc[1<<K];
void process_precalc(int i){for(size_t j=0;j<K;j++){int step=1<<j;
for(int k=0;k<step;k++){precalc[k+step]=precalc[k]^b[i+j];}}}
void solve(){int n,m,k;
cin>>n>>m>>k;
cp_algo::checkpoint("init");
string row;
for(int i=0;i<n;i++){cin>>row;
a[i]=row;}
for(auto&it:b){it=bitpack(k);}
for(int i=0;i<m;i++){cin>>row;
b[i]=row;}
for(auto&it:c){it=bitpack(k);}
for(auto&it:precalc){it=bitpack(k);}
cp_algo::checkpoint("read");
const int width=bitpack::width;
for(int j=0;j<m;j+=width){for(int offset=0;offset<width;offset+=K){process_precalc(j+offset);
for(int i=0;i<n;i++){c[i]^=precalc[uint8_t(a[i].word(j/width)>>offset)];}}}
cp_algo::checkpoint("mul");
for(int i=0;i<n;i++){cout<<c[i].to_string()<<"\n";}
cp_algo::checkpoint("write");
cp_algo::checkpoint<1>();}
signed main(){ios::sync_with_stdio(0);
cin.tie(0);
int t=1;
while(t--){solve();}}