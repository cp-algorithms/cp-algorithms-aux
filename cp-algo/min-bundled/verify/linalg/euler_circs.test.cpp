#line 1 "verify/linalg/euler_circs.test.cpp"
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
#line 1 "cp-algo/math/combinatorics.hpp"
#line 1 "cp-algo/math/common.hpp"
#line 6 "cp-algo/math/common.hpp"
namespace cp_algo::math{#ifdef CP_ALGO_MAXN
const int maxn=CP_ALGO_MAXN;
#else
const int maxn=1<<19;
#endif
const int magic=64;
auto bpow(auto const&x,auto n,auto const&one,auto op){if(n==0){return one;}else{auto t=bpow(x,n/2,one,op);
t=op(t,t);
if(n%2){t=op(t,x);}
return t;}}
auto bpow(auto x,auto n,auto ans){return bpow(x,n,ans,std::multiplies{});}
template<typename T>
T bpow(T const&x,auto n){return bpow(x,n,T(1));}
inline constexpr auto inv2(auto x){assert(x%2);
std::make_unsigned_t<decltype(x)>y=1;
while(y*x!=1){y*=2-x*y;}
return y;}}
#line 5 "cp-algo/math/combinatorics.hpp"
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
#line 1 "cp-algo/linalg/matrix.hpp"
#line 1 "cp-algo/random/rng.hpp"
#line 5 "cp-algo/random/rng.hpp"
namespace cp_algo::random{std::mt19937_64 gen(
std::chrono::steady_clock::now().time_since_epoch().count()
);
uint64_t rng(){return gen();}}
#line 1 "cp-algo/linalg/vector.hpp"
#line 1 "cp-algo/number_theory/modint.hpp"
#line 6 "cp-algo/number_theory/modint.hpp"
namespace cp_algo::math{template<typename modint,typename _Int>
struct modint_base{using Int=_Int;
using UInt=std::make_unsigned_t<Int>;
static constexpr size_t bits=sizeof(Int)*8;
using Int2=std::conditional_t<bits<=32,int64_t,__int128_t>;
using UInt2=std::conditional_t<bits<=32,uint64_t,__uint128_t>;
constexpr static Int mod(){return modint::mod();}
constexpr static Int remod(){return modint::remod();}
constexpr static UInt2 modmod(){return UInt2(mod())*mod();}
constexpr modint_base()=default;
constexpr modint_base(Int2 rr){to_modint().setr(UInt((rr+modmod())%mod()));}
modint inv()const{return bpow(to_modint(),mod()-2);}
modint operator-()const{modint neg;
neg.r=std::min(-r,remod()-r);
return neg;}
modint&operator/=(const modint&t){return to_modint()*=t.inv();}
modint&operator*=(const modint&t){r=UInt(UInt2(r)*t.r%mod());
return to_modint();}
modint&operator+=(const modint&t){r+=t.r;r=std::min(r,r-remod());
return to_modint();}
modint&operator-=(const modint&t){r-=t.r;r=std::min(r,r+remod());
return to_modint();}
modint operator+(const modint&t)const{return modint(to_modint())+=t;}
modint operator-(const modint&t)const{return modint(to_modint())-=t;}
modint operator*(const modint&t)const{return modint(to_modint())*=t;}
modint operator/(const modint&t)const{return modint(to_modint())/=t;}
auto operator==(const modint&t)const{return to_modint().getr()==t.getr();}
auto operator!=(const modint&t)const{return to_modint().getr()!=t.getr();}
auto operator<=(const modint&t)const{return to_modint().getr()<=t.getr();}
auto operator>=(const modint&t)const{return to_modint().getr()>=t.getr();}
auto operator<(const modint&t)const{return to_modint().getr()<t.getr();}
auto operator>(const modint&t)const{return to_modint().getr()>t.getr();}
Int rem()const{UInt R=to_modint().getr();
return R-(R>(UInt)mod()/2)*mod();}
constexpr void setr(UInt rr){r=rr;}
constexpr UInt getr()const{return r;}
static UInt modmod8(){return UInt(8*modmod());}
void add_unsafe(UInt t){r+=t;}
void pseudonormalize(){r=std::min(r,r-modmod8());}
modint const&normalize(){if(r>=(UInt)mod()){r%=mod();}
return to_modint();}
void setr_direct(UInt rr){r=rr;}
UInt getr_direct()const{return r;}
protected:UInt r;
private:constexpr modint&to_modint(){return static_cast<modint&>(*this);}
constexpr modint const&to_modint()const{return static_cast<modint const&>(*this);}};
template<typename modint>
concept modint_type=std::is_base_of_v<modint_base<modint,typename modint::Int>,modint>;
template<modint_type modint>
decltype(std::cin)&operator>>(decltype(std::cin)&in,modint&x){typename modint::UInt r;
auto&res=in>>r;
x.setr(r);
return res;}
template<modint_type modint>
decltype(std::cout)&operator<<(decltype(std::cout)&out,modint const&x){return out<<x.getr();}
template<auto m>
struct modint:modint_base<modint<m>,decltype(m)>{using Base=modint_base<modint<m>,decltype(m)>;
using Base::Base;
static constexpr Base::Int mod(){return m;}
static constexpr Base::UInt remod(){return m;}
auto getr()const{return Base::r;}};
template<typename Int=int>
struct dynamic_modint:modint_base<dynamic_modint<Int>,Int>{using Base=modint_base<dynamic_modint<Int>,Int>;
using Base::Base;
static Base::UInt m_reduce(Base::UInt2 ab){if(mod()%2==0)[[unlikely]]{return typename Base::UInt(ab%mod());}else{typename Base::UInt2 m=typename Base::UInt(ab)*imod();
return typename Base::UInt((ab+m*mod())>>Base::bits);}}
static Base::UInt m_transform(Base::UInt a){if(mod()%2==0)[[unlikely]]{return a;}else{return m_reduce(a*pw128());}}
dynamic_modint&operator*=(const dynamic_modint&t){Base::r=m_reduce(typename Base::UInt2(Base::r)*t.r);
return*this;}
void setr(Base::UInt rr){Base::r=m_transform(rr);}
Base::UInt getr()const{typename Base::UInt res=m_reduce(Base::r);
return std::min(res,res-mod());}
static Int mod(){return m;}
static Int remod(){return 2*m;}
static Base::UInt imod(){return im;}
static Base::UInt2 pw128(){return r2;}
static void switch_mod(Int nm){m=nm;
im=m%2?inv2(-m):0;
r2=static_cast<Base::UInt>(static_cast<Base::UInt2>(-1)%m+1);}
auto static with_mod(Int tmp,auto callback){struct scoped{Int prev=mod();
~scoped(){switch_mod(prev);}}_;
switch_mod(tmp);
return callback();}
private:static thread_local Int m;
static thread_local Base::UInt im,r2;};
template<typename Int>
Int thread_local dynamic_modint<Int>::m=1;
template<typename Int>
dynamic_modint<Int>::Base::UInt thread_local dynamic_modint<Int>::im=-1;
template<typename Int>
dynamic_modint<Int>::Base::UInt thread_local dynamic_modint<Int>::r2=0;}
#line 1 "cp-algo/util/big_alloc.hpp"
#line 7 "cp-algo/util/big_alloc.hpp"
#if defined(__linux__) || defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#  define CP_ALGO_USE_MMAP 1
#  include <sys/mman.h>
#else
#  define CP_ALGO_USE_MMAP 0
#endif
namespace cp_algo{template<typename T,std::size_t Align=32>
class big_alloc{static_assert(Align>=alignof(void*),"Align must be at least pointer-size");
static_assert(std::popcount(Align)==1,"Align must be a power of two");
public:using value_type=T;
template<class U>struct rebind{using other=big_alloc<U,Align>;};
constexpr bool operator==(const big_alloc&)const=default;
constexpr bool operator!=(const big_alloc&)const=default;
big_alloc()noexcept=default;
template<typename U,std::size_t A>
big_alloc(const big_alloc<U,A>&)noexcept{}
[[nodiscard]]T*allocate(std::size_t n){std::size_t padded=round_up(n*sizeof(T));
std::size_t align=std::max<std::size_t>(alignof(T),Align);
#if CP_ALGO_USE_MMAP
if(padded>=MEGABYTE){void*raw=mmap(nullptr,padded,
PROT_READ|PROT_WRITE,
MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
madvise(raw,padded,MADV_HUGEPAGE);
madvise(raw,padded,MADV_POPULATE_WRITE);
return static_cast<T*>(raw);}
#endif
return static_cast<T*>(::operator new(padded,std::align_val_t(align)));}
void deallocate(T*p,std::size_t n)noexcept{if(!p)return;
std::size_t padded=round_up(n*sizeof(T));
std::size_t align=std::max<std::size_t>(alignof(T),Align);
#if CP_ALGO_USE_MMAP
if(padded>=MEGABYTE){munmap(p,padded);return;}
#endif::operator delete(p,padded,std::align_val_t(align));}
private:static constexpr std::size_t MEGABYTE=1<<20;
static constexpr std::size_t round_up(std::size_t x)noexcept{return(x+Align-1)/Align*Align;}};
template<typename T>
using big_vector=std::vector<T,big_alloc<T>>;}
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
#line 1 "cp-algo/util/checkpoint.hpp"
#line 7 "cp-algo/util/checkpoint.hpp"
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
#line 15 "cp-algo/linalg/vector.hpp"
namespace cp_algo::linalg{template<typename base,class Alloc=big_alloc<base>>
struct vec:std::basic_string<base,std::char_traits<base>,Alloc>{using Base=std::basic_string<base,std::char_traits<base>,Alloc>;
using Base::Base;
vec(Base const&t):Base(t){}
vec(Base&&t):Base(std::move(t)){}
vec(size_t n):Base(n,base()){}
vec(auto&&r):Base(std::ranges::to<Base>(r)){}
static vec ei(size_t n,size_t i){vec res(n);
res[i]=1;
return res;}
auto operator-()const{return*this|std::views::transform([](auto x){return-x;});}
auto operator*(base t)const{return*this|std::views::transform([t](auto x){return x*t;});}
auto operator*=(base t){for(auto&it:*this){it*=t;}
return*this;}
virtual void add_scaled(vec const&b,base scale,size_t i=0){if(scale!=base(0)){for(;i<size(*this);i++){(*this)[i]+=scale*b[i];}}}
virtual vec const&normalize(){return static_cast<vec&>(*this);}
virtual base normalize(size_t i){return(*this)[i];}
void read(){for(auto&it:*this){std::cin>>it;}}
void print()const{for(auto&it:*this){std::cout<<it<<" ";}
std::cout<<"\n";}
static vec random(size_t n){vec res(n);
std::ranges::generate(res,random::rng);
return res;}
vec operator|(vec const&t)const{return std::views::join(std::array{std::views::all(*this),
std::views::all(t)});}
std::pair<size_t,base>find_pivot(){if(pivot==size_t(-1)){pivot=0;
while(pivot<size(*this)&&normalize(pivot)==base(0)){pivot++;}
if(pivot<size(*this)){pivot_inv=base(1)/(*this)[pivot];}}
return{pivot,pivot_inv};}
void reduce_by(vec&t){auto[pivot,pinv]=t.find_pivot();
if(pivot<size(*this)){add_scaled(t,-normalize(pivot)*pinv,pivot);}}
private:size_t pivot=-1;
base pivot_inv;};
template<math::modint_type base,class Alloc=big_alloc<base>>
struct modint_vec:vec<base,Alloc>{using Base=vec<base,Alloc>;
using Base::Base;
modint_vec(Base const&t):Base(t){}
modint_vec(Base&&t):Base(std::move(t)){}
void add_scaled(Base const&b,base scale,size_t i=0)override{static_assert(base::bits>=64,"Only wide modint types for linalg");
if(scale!=base(0)){assert(Base::size()==b.size());
size_t n=size(*this);
u64x4 scaler=u64x4()+scale.getr();
if(is_aligned(&(*this)[0])&&is_aligned(&b[0]))
for(i-=i%4;i+3<n;i+=4){auto&ai=vector_cast<u64x4>((*this)[i]);
auto bi=vector_cast<u64x4 const>(b[i]);
#ifdef __AVX2__
ai+=u64x4(_mm256_mul_epu32(__m256i(scaler),__m256i(bi)));
#else
ai+=scaler*bi;
#endif}
for(;i<n;i++){(*this)[i].add_unsafe(b[i].getr_direct()*scale.getr());}
if(++counter==4){for(auto&it:*this){it.pseudonormalize();}
counter=0;}}}
Base const&normalize()override{for(auto&it:*this){it.normalize();}
return*this;}
base normalize(size_t i)override{return(*this)[i].normalize();}
private:size_t counter=0;};}
#line 7 "cp-algo/linalg/matrix.hpp"
#include <optional>
#line 11 "cp-algo/linalg/matrix.hpp"
namespace cp_algo::linalg{enum gauss_mode{normal,reverse};
template<typename base_t,class _vec_t=std::conditional_t<
math::modint_type<base_t>,
modint_vec<base_t>,
vec<base_t>>>
struct matrix:std::vector<_vec_t>{using vec_t=_vec_t;
using base=base_t;
using Base=std::vector<vec_t>;
using Base::Base;
matrix(size_t n):Base(n,vec_t(n)){}
matrix(size_t n,size_t m):Base(n,vec_t(m)){}
matrix(Base const&t):Base(t){}
matrix(Base&&t):Base(std::move(t)){}
template<std::ranges::input_range R>
matrix(R&&r):Base(std::ranges::to<Base>(std::forward<R>(r))){}
size_t n()const{return size(*this);}
size_t m()const{return n()?size(row(0)):0;}
void resize(size_t n,size_t m){Base::resize(n);
for(auto&it:*this){it.resize(m);}}
auto&row(size_t i){return(*this)[i];}
auto const&row(size_t i)const{return(*this)[i];}
auto elements(){return*this|std::views::join;}
auto elements()const{return*this|std::views::join;}
matrix operator-()const{return*this|std::views::transform([](auto x){return vec_t(-x);});}
matrix&operator+=(matrix const&t){for(auto[a,b]:std::views::zip(elements(),t.elements())){a+=b;}
return*this;}
matrix&operator-=(matrix const&t){for(auto[a,b]:std::views::zip(elements(),t.elements())){a-=b;}
return*this;}
matrix operator+(matrix const&t)const{return matrix(*this)+=t;}
matrix operator-(matrix const&t)const{return matrix(*this)-=t;}
matrix&operator*=(base t){for(auto&it:*this)it*=t;return*this;}
matrix operator*(base t)const{return matrix(*this)*=t;}
matrix&operator/=(base t){return*this*=base(1)/t;}
matrix operator/(base t)const{return matrix(*this)/=t;}
matrix&operator*=(matrix const&t){return*this=*this*t;}
void read_transposed(){for(size_t j=0;j<m();j++){for(size_t i=0;i<n();i++){std::cin>>(*this)[i][j];}}}
void read(){for(auto&it:*this){it.read();}}
void print()const{for(auto const&it:*this){it.print();}}
static matrix block_diagonal(std::vector<matrix>const&blocks){size_t n=0;
for(auto&it:blocks){assert(it.n()==it.m());
n+=it.n();}
matrix res(n);
n=0;
for(auto&it:blocks){for(size_t i=0;i<it.n();i++){std::ranges::copy(it[i],begin(res[n+i])+n);}
n+=it.n();}
return res;}
static matrix random(size_t n,size_t m){matrix res(n,m);
std::ranges::generate(res,std::bind(vec_t::random,m));
return res;}
static matrix random(size_t n){return random(n,n);}
static matrix eye(size_t n){matrix res(n);
for(size_t i=0;i<n;i++){res[i][i]=1;}
return res;}
matrix operator|(matrix const&b)const{assert(n()==b.n());
matrix res(n(),m()+b.m());
for(size_t i=0;i<n();i++){res[i]=row(i)|b[i];}
return res;}
void assign_submatrix(auto viewx,auto viewy,matrix const&t){for(auto[a,b]:std::views::zip(*this|viewx,t)){std::ranges::copy(b,begin(a|viewy));}}
auto submatrix(auto viewx,auto viewy)const{return*this|viewx|std::views::transform([viewy](auto const&y){return y|viewy;});}
matrix T()const{matrix res(m(),n());
for(size_t i=0;i<n();i++){for(size_t j=0;j<m();j++){res[j][i]=row(i)[j];}}
return res;}
matrix operator*(matrix const&b)const{assert(m()==b.n());
matrix res(n(),b.m());
for(size_t i=0;i<n();i++){for(size_t j=0;j<m();j++){res[i].add_scaled(b[j],row(i)[j]);}}
return res.normalize();}
vec_t apply(vec_t const&x)const{return(matrix(1,x)**this)[0];}
matrix pow(uint64_t k)const{assert(n()==m());
return bpow(*this,k,eye(n()));}
matrix&normalize(){for(auto&it:*this){it.normalize();}
return*this;}
template<gauss_mode mode=normal>
void eliminate(size_t i,size_t k){auto kinv=base(1)/row(i).normalize()[k];
for(size_t j=(mode==normal)*i;j<n();j++){if(j!=i){row(j).add_scaled(row(i),-row(j).normalize(k)*kinv);}}}
template<gauss_mode mode=normal>
void eliminate(size_t i){row(i).normalize();
for(size_t j=(mode==normal)*i;j<n();j++){if(j!=i){row(j).reduce_by(row(i));}}}
template<gauss_mode mode=normal>
matrix&gauss(){for(size_t i=0;i<n();i++){eliminate<mode>(i);}
return normalize();}
template<gauss_mode mode=normal>
auto echelonize(size_t lim){return gauss<mode>().sort_classify(lim);}
template<gauss_mode mode=normal>
auto echelonize(){return echelonize<mode>(m());}
size_t rank()const{if(n()>m()){return T().rank();}
return size(matrix(*this).echelonize()[0]);}
base det()const{assert(n()==m());
matrix b=*this;
b.echelonize();
base res=1;
for(size_t i=0;i<n();i++){res*=b[i][i];}
return res;}
std::pair<base,matrix>inv()const{assert(n()==m());
matrix b=*this|eye(n());
if(size(b.echelonize<reverse>(n())[0])<n()){return{0,{}};}
base det=1;
for(size_t i=0;i<n();i++){det*=b[i][i];
b[i]*=base(1)/b[i][i];}
return{det,b.submatrix(std::views::all,std::views::drop(n()))};}
auto kernel()const{auto A=*this;
auto[pivots,free]=A.template echelonize<reverse>();
matrix sols(size(free),m());
for(size_t j=0;j<size(pivots);j++){base scale=base(1)/A[j][pivots[j]];
for(size_t i=0;i<size(free);i++){sols[i][pivots[j]]=A[j][free[i]]*scale;}}
for(size_t i=0;i<size(free);i++){sols[i][free[i]]=-1;}
return sols;}
std::optional<std::array<matrix,2>>solve(matrix t)const{matrix sols=(*this|t).kernel();
if(sols.n()<t.m()||matrix(sols.submatrix(
std::views::drop(sols.n()-t.m()),
std::views::drop(m())
))!=-eye(t.m())){return std::nullopt;}else{return std::array{matrix(sols.submatrix(std::views::drop(sols.n()-t.m()),std::views::take(m()))),
matrix(sols.submatrix(std::views::take(sols.n()-t.m()),std::views::take(m())))};}}
auto sort_classify(size_t lim){size_t rk=0;
std::vector<size_t>free,pivots;
for(size_t j=0;j<lim;j++){for(size_t i=rk+1;i<n()&&row(rk)[j]==base(0);i++){if(row(i)[j]!=base(0)){std::swap(row(i),row(rk));
row(rk)=-row(rk);}}
if(rk<n()&&row(rk)[j]!=base(0)){pivots.push_back(j);
rk++;}else{free.push_back(j);}}
return std::array{pivots,free};}};
template<typename base_t>
auto operator*(base_t t,matrix<base_t>const&A){return A*t;}}
#line 8 "verify/linalg/euler_circs.test.cpp"
using namespace std;
using namespace cp_algo::math;
using namespace cp_algo::linalg;
const int64_t mod=998244353;
using base=modint<mod>;
void solve(){int n,m;
cin>>n>>m;
matrix<base>a(n);
int r=0;
vector<int>indeg(n),outdeg(n);
for(int i=0;i<m;i++){int u,v;
cin>>u>>v;
a[u][v]-=1;
a[v][v]+=1;
outdeg[u]++;
indeg[v]++;
r=v;}
a[r][r]=fact<base>(indeg[r]-1);
for(int i=0;i<n;i++){if(i==r){continue;}
if(indeg[i]!=outdeg[i]){cout<<0<<"\n";
return;}
if(indeg[i]!=0){a[i]*=fact<base>(indeg[i]-1);}else{a[i][i]=1;}
a[r][i]=a[i][r]=0;}
cout<<a.det()<<"\n";}
signed main(){ios::sync_with_stdio(0);
cin.tie(0);
int t=1;
while(t--){solve();}}