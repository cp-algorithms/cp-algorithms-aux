#line 1 "verify/structures/fenwick/ordered_set.test.cpp"
#define PROBLEM "https:#pragma GCC optimize("Ofast,unroll-loops")
#include <bits/stdc++.h>
#line 1 "blazingio/blazingio.min.hpp"
#define M$(x,...)_mm256_##x##_epi8(__VA_ARGS__)
#define $u(...)__VA_ARGS__
#if __APPLE__
#define $m(A,B)A
#else
#define $m(A,B)B
#endif
#if _WIN32
#define $w(A,B)A
#else
#define $w(A,B)B
#endif
#if __i386__|_M_IX86
#define $H(A,B)A
#else
#define $H(A,B)B
#endif
#if __aarch64__
#define $a(A,B)A
#else
#define $a(A,B)B
#endif
#define $P(x)void F(x K){#define $T template<$c T
#define $c class
#define $C constexpr
#define $R return
#define $O operator
#define u$ uint64_t
#define $r $R*this;
#line 41 "blazingio/blazingio.min.hpp"
#include $a(<arm_neon.h>,<immintrin.h>)
#line 43 "blazingio/blazingio.min.hpp"
#include $w(<windows.h>,<sys/mman.h>)
#include<sys/stat.h>
#include $w(<io.h>,<unistd.h>)
#include $w(<ios>,<sys/resource.h>)
#if _MSC_VER
#define __builtin_add_overflow(a,b,c)_addcarry_u64(0,a,b,c)
#define $s
#else
$H(,u$ _umul128(u$ a,u$ b,u$*D){auto x=(__uint128_t)a*b;*D=u$(x>>64);$R(u$)x;})
#define $s $a(,__attribute__((target("avx2"))))
#endif
#define $z $a(16,32)
#define $t $a(uint8x16_t,__m256i)
#define $I $w(__forceinline,__attribute__((always_inline)))
#define $F M(),
#define E$(x)if(!(x))abort();
$w(LONG WINAPI $x(_EXCEPTION_POINTERS*);,)namespace $f{using namespace std;struct B{enum $c A:char{}c;B&$O=(char x){c=A{x};$r}$O char(){$R(char)c;}};$C u$ C=~0ULL/255;struct D{string&K;};static B E[65568];template<int F>struct G{B*H,*S;void K(off_t C){$w(char*D=(char*)VirtualAlloc(0,(C+8191)&-4096,8192,1);E$(D)E$(VirtualFree(D,0,32768))DWORD A=C&-65536;E$(!A||MapViewOfFileEx(CreateFileMapping(GetStdHandle(-10),0,2,0,A,0),4,0,0,0,D)==D)E$(VirtualAlloc(D+A,65536,12288,4)==D+A)E$(~_lseek(0,A,0))DWORD E=0;ReadFile(GetStdHandle(-10),D+A,65536,&E,0);,int A=getpagesize();char*D=(char*)mmap(0,C+A,3,2,0,0);E$(D!=(void*)-1)E$(mmap(D+((C+A-1)&-A),A,3,$m(4114,50),-1,0)!=(void*)-1))H=(B*)D+C;*H=10;H[1]=48;H[2]=0;S=(B*)D;}void L(){H=S=E;}$I void M(){if(F&&S==H){$w(DWORD A=0;ReadFile(GetStdHandle(-10),S=E,65536,&A,0);,$a($u(register long A asm("x0")=0,D asm("x1")=(long)E,G asm("x2")=65536,C asm($m("x16","x8"))=$m(3,63);asm volatile("svc 0"$m("x80",):"+r"(A),"+r"(D):"r"(C),"r"(G));S=launder(E);),off_t A=$H(3,$m(33554435,0));B*D=E;asm volatile($H("int $128","syscall"):"+a"(A),$H("+c"(D):"b","+S"(D):"D")(0),"d"(65536)$H(,$u(:"rcx","r11")));S=D;))H=S+A;*H=10;if(!A)E[1]=48,E[2]=0;}}$T>$I void N(T&x){while($F(*S&240)==48)x=T(x*10+(*S++-48));}$T>$I decltype((void)~T{1})O(T&x){M();int A=is_signed_v<T>&&*S==45;S+=A;N(x=0);x=A?1+~x:x;}$T>$I decltype((void)T{1.})O(T&x){M();int A=*S==45;S+=A;$F S+=*S==43;u$ n=0;int i=0;for(;i<18&&($F*S&240)==48;i++)n=n*10+*S++-48;int B=20;int C=*S==46;S+=C;for(;i<18&&($F*S&240)==48;i++)n=n*10+*S++-48,B-=C;x=(T)n;while(($F*S&240)==48)x=x*10+*S++-48,B-=C;if(*S==46)S++,C=1;while(($F*S&240)==48)x=x*10+*S++-48,B-=C;int D;if((*S|32)==101)S++,$F S+=*S==43,O(D),B+=D;static $C auto E=[](){array<T,41>E{};T x=1;for(int i=21;i--;)E[40-i]=x,E[i]=1/x,x*=10;$R E;}();while(B>40)x*=(T)1e10,B-=10;while(B<0)x*=(T)1e-10,B+=10;x*=E[B];x=A?-x:x;}$I void O(bool&x){$F x=*S++==49;}$I void O(char&x){$F x=*S++;}$I void O(uint8_t&x){$F x=*S++;}$I void O(int8_t&x){$F x=*S++;}$T>$s void P(string&K,T C){M();B*G=S;C();K.assign((char*)G,S-G);while(F&&S==H&&($F H!=E)){C();K.append(E,S);}}$s void O(string&K){P(K,[&]()$s{B*p=S;$w(ULONG R;,)$t x;$a(uint64x2_t A;while(memcpy(&x,p,16),A=uint64x2_t(x<33),!(A[0]|A[1]))p+=16;S=p+(A[0]?0:8)+$w((_BitScanForward64(&R,A[0]?A[0]:A[1]),R),__builtin_ctzll(A[0]?A[0]:A[1]))/8;,int J;$t C=M$(set1,32);while(memcpy(&x,p,32),!(J=M$(movemask,M$(cmpeq,C,_mm256_max_epu8(C,x)))))p+=32;S=p+$w((_BitScanForward(&R,J),R),__builtin_ctz(J));)});}$s void O(D&A){P(A.K,[&](){S=(B*)memchr(S,10,H-S+1);});if(A.K.size()&&A.K.back()==13)A.K.pop_back();if(A.K.empty()||S<H)S+=*S==10;}$T>$I void O(complex<T>&K){T A,B{};if($F*S==40){S++;O(A);if($F*S++==44)Q(B),S++;}else O(A);K={A,B};}template<size_t N>$s void O(bitset<N>&K){if(N>4095&&!*this)$R;ptrdiff_t i=N;while(i)if($F i%$z||H-S<$z)K[--i]=*S++==49;else{B*p=S;for(int64_t j=0;j<min(i,H-S)/$z;j++){i-=$z;$t x;memcpy(&x,p,$z);$a(auto B=(uint8x16_t)vdupq_n_u64(~2ULL/254)&(48-x);auto C=vzip_u8(vget_high_u8(B),vget_low_u8(B));auto y=vaddvq_u16((uint16x8_t)vcombine_u8(C.val[0],C.val[1]));,u$ a=~0ULL/65025;auto y=$w(_byteswap_ulong,__builtin_bswap32)(M$(movemask,M$(shuffle,_mm256_slli_epi32(x,7),_mm256_set_epi64x(a+C*24,a+C*16,a+C*8,a))));)p+=$z;memcpy((char*)&K+i/8,&y,$z/8);}S=p;}}$T>$I void Q(T&K){if(!is_same_v<T,D>)while($F(uint8_t)*S<33)S++;O(K);}$O bool(){$R!!*this;}bool $O!(){$R S>H;}};struct U{G<0>A;G<1>B;U(){struct stat D;E$(~fstat(0,&D))(D.st_mode>>12)==8?A.K(D.st_size):B.L();}U*tie(nullptr_t){$R this;}void sync_with_stdio(bool){}$T>$I U&$O>>(T&K){A.S?A.Q(K):B.Q(K);$r}$O bool(){$R!!*this;}bool $O!(){$R A.S?!A:!B;}};short A[100];char L[64]{1};struct
V{char*D;B*S;int J;V(){$w(E$(D=(char*)VirtualAlloc(0,536870912,8192,4))E$(VirtualAlloc(D,4096,4096,260))AddVectoredExceptionHandler(1,$x);,size_t C=536870912;$m(,rlimit E;getrlimit(RLIMIT_AS,&E);if(~E.rlim_cur)C=25165824;)D=(char*)mmap(0,C,3,$m(4162,16418),-1,0);E$(D!=(void*)-1))S=(B*)D;for(int i=0;i<100;i++)A[i]=short((48+i/10)|((48+i%10)<<8));for(int i=1;i<64;i++)L[i]=L[i-1]+(0x8922489224892249>>i&1);}~V(){flush($w(!J,));}void flush($w(int F=0,)){$w(J=1;auto E=GetStdHandle(-11);auto C=F?ReOpenFile(E,1073741824,7,2684354560):(void*)-1;DWORD A;E$(C==(void*)-1?WriteFile(E,D,DWORD((char*)S-D),&A,0):(WriteFile(C,D,DWORD(((char*)S-D+4095)&-4096),&A,0)&&~_chsize(1,int((char*)S-D)))),auto G=D;ssize_t A;while((A=write(1,G,(char*)S-G))>0)G+=A;E$(~A))S=(B*)D;}$P(char)*S++=K;}$P(uint8_t)*S++=K;}$P(int8_t)*S++=K;}$P(bool)*S++=48+K;}$T>decltype((void)~T{1})F(T K){using D=make_unsigned_t<T>;D C=K;if(K<0)F('-'),C=1+~C;static $C auto N=[](){array<D,5*sizeof(T)/2>N{};D n=1;for(size_t i=1;i<N.size();i++)n*=10,N[i]=n;$R N;}();$w(ULONG M;,)int G=L[$w(($H(_BitScanReverse(&M,ULONG((int64_t)C>>32))?M+=32:_BitScanReverse(&M,(ULONG)C|1),_BitScanReverse64(&M,C|1)),M),63^__builtin_clzll(C|1))];G-=C<N[G-1];short H[20];if $C(sizeof(T)==2){auto n=33555U*C-C/2;u$ H=A[n>>25];n=(n&33554431)*25;H|=A[n>>23]<<16;H|=u$(48+((n&8388607)*5>>22))<<32;H>>=40-G*8;memcpy(S,&H,8);}else if $C(sizeof(T)==4){auto n=1441151881ULL*C;$H(n>>=25;n++;for(int i=0;i<5;i++){H[i]=A[n>>32];n=(n&~0U)*100;},int K=57;auto J=~0ULL>>7;for(int i=0;i<5;i++){H[i]=A[n>>K];n=(n&J)*25;K-=2;J/=4;})memcpy(S,(B*)H+10-G,16);}else{$H($u(if(C<(1ULL<<32)){$R F((uint32_t)C);}auto J=(u$)1e10;auto x=C/J,y=C%J;int K=100000,b[]{int(x/K),int(x%K),int(y/K),int(y%K)};B H[40];for(int i=0;i<4;i++){int n=int((429497ULL*b[i]>>7)+1);B*p=H+i*5;*p=48+char(n>>25);n=(n&~0U>>7)*25;memcpy(p+1,A+(n>>23),2);memcpy(p+3,A+((n&~0U>>9)*25>>21),2);}),$u(u$ D,E=_umul128(18,C,&D),F;_umul128(0x725dd1d243aba0e8,C,&F);D+=__builtin_add_overflow(E,F+1,&E);for(int i=0;i<10;i++)H[i]=A[D],E=_umul128(100,E,&D);))memcpy(S,(B*)H+20-G,20);}S+=G;}$T>decltype((void)T{1.})F(T K){if(K<0)F('-'),K=-K;auto G=[&](){auto x=u$(K*1e12);$H($u(x-=x>999999999999;uint32_t n[]{uint32_t(x/1000000*429497>>7)+1,uint32_t(x%1000000*429497>>7)+1};int K=25,J=~0U>>7;for(int i=0;i<3;i++){for(int j=0;j<2;j++)memcpy(S+i*2+j*6,A+(n[j]>>K),2),n[j]=(n[j]&J)*25;K-=2;J/=4;}S+=12;),$u(u$ D,E=_umul128(472236648287,x,&D)>>8;E|=D<<56;D>>=8;E++;for(int i=0;i<6;i++)memcpy(S,A+D,2),S+=2,E=_umul128(100,E,&D);))};if(K==0)$R F('0');if(K>=1e16){K*=(T)1e-16;int B=16;while(K>=1)K*=(T).1,B++;F("0.");G();F('e');F(B);}else if(K>=1){auto B=(u$)K;F(B);if((K-=(T)B)>0)F('.'),G();}else F("0."),G();}$P(const char*)$w(size_t A=strlen(K);memcpy((char*)S,K,A);S+=A;,S=(B*)stpcpy((char*)S,K);)}$P(const uint8_t*)F((char*)K);}$P(const int8_t*)F((char*)K);}$P(string_view)memcpy(S,K.data(),K.size());S+=K.size();}$T>$P(complex<T>)*this<<'('<<K.real()<<','<<K.imag()<<')';}template<size_t N>$s $P(const bitset<N>&)auto i=N;while(i%$z)*S++=48+K[--i];B*p=S;while(i){i-=$z;$a(short,int)x;memcpy(&x,(char*)&K+i/8,$z/8);$a(auto A=(uint8x8_t)vdup_n_u16(x);vst1q_u8((uint8_t*)p,48-vtstq_u8(vcombine_u8(vuzp2_u8(A,A),vuzp1_u8(A,A)),(uint8x16_t)vdupq_n_u64(~2ULL/254)));,auto b=_mm256_set1_epi64x(~2ULL/254);_mm256_storeu_si256(($t*)p,M$(sub,M$(set1,48),M$(cmpeq,_mm256_and_si256(M$(shuffle,_mm256_set1_epi32(x),_mm256_set_epi64x(0,C,C*2,C*3)),b),b)));)p+=$z;}S=p;}$T>V&$O<<(const T&K){F(K);$r}V&$O<<(V&(*A)(V&)){$R A(*this);}};struct W{$T>W&$O<<(const T&K){$r}W&$O<<(W&(*A)(W&)){$R A(*this);}};}namespace std{$f::U i$;$f::V o$;$f::W e$;$f::U&getline($f::U&B,string&K){$f::D A{K};$R B>>A;}$f::V&flush($f::V&B){if(!i$.A.S)B.flush();$R B;}$f::V&endl($f::V&B){$R B<<'\n'<<flush;}$f::W&endl($f::W&B){$R B;}$f::W&flush($f::W&B){$R B;}}$w(LONG WINAPI $x(_EXCEPTION_POINTERS*A){auto C=A->ExceptionRecord;auto B=C->ExceptionInformation[1];if(C->ExceptionCode==2147483649&&B-(ULONG_PTR)std::o$.D<0x40000000){E$(VirtualAlloc((char*)B,16777216,4096,4)&&VirtualAlloc((char*)(B+16777216),4096,4096,260))$R-1;}$R 0;},)
#define freopen(...)if(freopen(__VA_ARGS__)==stdin)std::i$=$f::U{}
#define cin i$
#define cout o$
#ifdef ONLINE_JUDGE
#define cerr e$
#define clog e$
#endif
#line 1 "cp-algo/structures/fenwick_set.hpp"
#line 1 "cp-algo/structures/fenwick.hpp"
#line 5 "cp-algo/structures/fenwick.hpp"
namespace cp_algo::structures{template<typename Op>
struct inverse_op{};
template<typename T>
struct inverse_op<std::plus<T>>{static T apply(T const&a,T const&b){return a-b;}};
template<typename T>
struct inverse_op<std::multiplies<T>>{static T apply(T const&a,T const&b){return a/b;}};
template<typename T,std::ranges::range Container=std::vector<T>,typename Op=std::plus<T>>
struct fenwick{Op op;
size_t n;
Container data;
fenwick(auto&&range,Op&&op=Op{}):op(std::move(op)){assign(std::move(range));}
void to_prefix_folds(){for(size_t i=1;i<n;i++){if(i+(i&-i)<=n){data[i+(i&-i)]=op(data[i+(i&-i)],data[i]);}}}
void assign(auto&&range){n=size(range)-1;
data=std::move(range);
to_prefix_folds();}
void update(size_t x,T const&v){for(++x;x<=n;x+=x&-x){data[x]=op(data[x],v);}}
T prefix_fold(size_t r)const{assert(r<=n);
T res={};
for(;r;r-=r&-r){res=op(res,data[r]);}
return res;}
T range_fold(size_t l,size_t r)const{return inverse_op<Op>::apply(prefix_fold(r),prefix_fold(l));}
auto prefix_lower_bound(T k)const{size_t x=0;
T pref={};
for(size_t i=std::bit_floor(n);i;i/=2){if(x+i<=n&&op(pref,data[x+i])<=k){pref=op(pref,data[x+i]);
x+=i;}}
return std::pair{x,pref};}};
template<std::ranges::range Container,typename Op>
fenwick(Container&&,Op&&)->fenwick<std::ranges::range_value_t<Container>,Container,Op>;
template<std::ranges::range Container>
fenwick(Container&&)->fenwick<std::ranges::range_value_t<Container>,Container>;
auto maxer=[](auto const&a,auto const&b){return std::max(a,b);};
template<typename T,std::ranges::range Container=std::vector<T>>
struct fenwick_max:fenwick<T,Container,decltype(maxer)>{using fenwick<T,Container,decltype(maxer)>::fenwick;};
template<std::ranges::range Container>
fenwick_max(Container&&)->fenwick_max<std::ranges::range_value_t<Container>,Container>;}
#line 1 "cp-algo/structures/bit_array.hpp"
#line 1 "cp-algo/util/bit.hpp"
#line 1 "cp-algo/util/simd.hpp"
#include <experimental/simd>
#line 7 "cp-algo/util/simd.hpp"
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
#line 6 "cp-algo/util/bit.hpp"
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
#line 5 "cp-algo/structures/bit_array.hpp"
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
#line 5 "cp-algo/structures/fenwick_set.hpp"
namespace cp_algo::structures{template<size_t maxc>
using popcount_array=std::array<int,maxc/bit_width<uint64_t>+1>;
template<size_t maxc>
struct fenwick_set:fenwick<int,popcount_array<maxc>>{using Base=fenwick<int,popcount_array<maxc>>;
static constexpr size_t word=bit_width<uint64_t>;
size_t sz=0;
bit_array<maxc>bits;
fenwick_set():Base(popcount_array<maxc>{}){}
fenwick_set(auto&&range):fenwick_set(){for(auto x:range){Base::data[x/word+1]+=1;
if(!bits.test(x)){sz++;
bits.flip(x);}}
Base::to_prefix_folds();}
void insert(size_t x){if(bits.test(x))return;
Base::update(x/word,1);
bits.flip(x);
sz++;}
void erase(size_t x){if(!bits.test(x))return;
Base::update(x/word,-1);
bits.flip(x);
sz--;}
size_t order_of_key(size_t x)const{return Base::prefix_fold(x/word)+order_of_bit(bits.word(x/word),x%word);}
size_t find_by_order(size_t order)const{if(order>=sz){return-1;}
auto[x,pref]=Base::prefix_lower_bound((int)order);
return x*word+kth_set_bit(bits.word(x),order-pref);}
size_t lower_bound(size_t x)const{if(bits.test(x)){return x;}
auto order=order_of_key(x);
return order<sz?find_by_order(order):-1;}
size_t pre_upper_bound(size_t x)const{if(bits.test(x)){return x;}
auto order=order_of_key(x);
return order?find_by_order(order-1):-1;}};}
#line 1 "cp-algo/util/compress_coords.hpp"
#line 1 "cp-algo/util/sort.hpp"
#line 6 "cp-algo/util/sort.hpp"
#include <ranges>
#line 8 "cp-algo/util/sort.hpp"
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
#line 5 "cp-algo/util/compress_coords.hpp"
namespace cp_algo{auto compress_coords(auto&&coords){using T=std::decay_t<std::unwrap_reference_t<
std::ranges::range_value_t<decltype(coords)>
>>;
std::vector<T>original;
if(empty(coords)){return original;}
original.reserve(size(coords));
radix_sort(coords);
int idx=-1;
T prev=~coords.front();
for(auto&x:coords){if(x!=prev){idx++;
prev=x;
original.push_back(x);}
x.get()=idx;}
return original;}}
#line 8 "verify/structures/fenwick/ordered_set.test.cpp"
using namespace std;
using cp_algo::structures::fenwick_set;
void solve(){int n,q;
cin>>n>>q;
vector a(n,0);
vector<reference_wrapper<int>>coords;
for(auto&it:a){cin>>it;
coords.push_back(ref(it));}
vector queries(q,pair{0,0});
for(auto&[t,x]:queries){cin>>t>>x;
if(t!=2){coords.push_back(ref(x));}}
auto values=cp_algo::compress_coords(coords);
const int maxc=1e6;
fenwick_set<maxc>me(a);
for(auto[t,x]:queries){if(t==0){me.insert(x);}else if(t==1){me.erase(x);}else if(t==2){auto res=(int)me.find_by_order(x-1);
cout<<(res==-1?-1:values[res])<<'\n';}else if(t==3){cout<<me.order_of_key(x+1)<<'\n';}else if(t==4){auto res=(int)me.pre_upper_bound(x);
cout<<(res==-1?-1:values[res])<<'\n';}else if(t==5){auto res=(int)me.lower_bound(x);
cout<<(res==-1?-1:values[res])<<'\n';}}}
signed main(){ios::sync_with_stdio(0);
cin.tie(0);
int t=1;
while(t--){solve();}}