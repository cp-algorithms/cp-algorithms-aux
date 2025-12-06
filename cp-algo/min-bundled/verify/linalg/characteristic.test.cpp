#line 1 "verify/linalg/characteristic.test.cpp"
#define PROBLEM "https:#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_MAXN 1 << 10
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
#line 1 "cp-algo/linalg/frobenius.hpp"
#line 1 "cp-algo/math/poly.hpp"
#line 1 "cp-algo/math/poly/impl/euclid.hpp"
#line 1 "cp-algo/math/affine.hpp"
#include <optional>
#line 7 "cp-algo/math/affine.hpp"
namespace cp_algo::math{template<typename base>
struct lin{base a=1,b=0;
std::optional<base>c;
lin(){}
lin(base b):a(0),b(b){}
lin(base a,base b):a(a),b(b){}
lin(base a,base b,base _c):a(a),b(b),c(_c){}
lin operator*(const lin&t){assert(c&&t.c&&*c==*t.c);
return{a*t.b+b*t.a,b*t.b+a*t.a*(*c),*c};}
lin apply(lin const&t)const{return{a*t.a,a*t.b+b};}
void prepend(lin const&t){*this=t.apply(*this);}
base eval(base x)const{return a*x+b;}};
template<typename base>
struct linfrac{base a,b,c,d;
linfrac():a(1),b(0),c(0),d(1){}
linfrac(base a):a(a),b(1),c(1),d(0){}
linfrac(base a,base b,base c,base d):a(a),b(b),c(c),d(d){}
linfrac operator*(linfrac t)const{return t.prepend(linfrac(*this));}
linfrac operator-()const{return{-a,-b,-c,-d};}
linfrac adj()const{return{d,-b,-c,a};}
linfrac&prepend(linfrac const&t){t.apply(a,c);
t.apply(b,d);
return*this;}
void apply(base&A,base&B)const{std::tie(A,B)=std::pair{a*A+b*B,c*A+d*B};}};}
#line 1 "cp-algo/math/fft.hpp"
#line 1 "cp-algo/number_theory/modint.hpp"
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
#line 1 "cp-algo/random/rng.hpp"
#line 5 "cp-algo/random/rng.hpp"
namespace cp_algo::random{std::mt19937_64 gen(
std::chrono::steady_clock::now().time_since_epoch().count()
);
uint64_t rng(){return gen();}}
#line 1 "cp-algo/math/cvector.hpp"
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
#line 1 "cp-algo/util/complex.hpp"
#line 5 "cp-algo/util/complex.hpp"
namespace cp_algo{template<typename T>
struct complex{using value_type=T;
T x,y;
constexpr complex():x(),y(){}
constexpr complex(T x):x(x),y(){}
constexpr complex(T x,T y):x(x),y(y){}
complex&operator*=(T t){x*=t;y*=t;return*this;}
complex&operator/=(T t){x/=t;y/=t;return*this;}
complex operator*(T t)const{return complex(*this)*=t;}
complex operator/(T t)const{return complex(*this)/=t;}
complex&operator+=(complex t){x+=t.x;y+=t.y;return*this;}
complex&operator-=(complex t){x-=t.x;y-=t.y;return*this;}
complex operator*(complex t)const{return{x*t.x-y*t.y,x*t.y+y*t.x};}
complex operator/(complex t)const{return*this*t.conj()/t.norm();}
complex operator+(complex t)const{return complex(*this)+=t;}
complex operator-(complex t)const{return complex(*this)-=t;}
complex&operator*=(complex t){return*this=*this*t;}
complex&operator/=(complex t){return*this=*this/t;}
complex operator-()const{return{-x,-y};}
complex conj()const{return{x,-y};}
T norm()const{return x*x+y*y;}
T abs()const{return std::sqrt(norm());}
T const real()const{return x;}
T const imag()const{return y;}
T&real(){return x;}
T&imag(){return y;}
static constexpr complex polar(T r,T theta){return{T(r*cos(theta)),T(r*sin(theta))};}
auto operator<=>(complex const&t)const=default;};
template<typename T>
complex<T>operator*(auto x,complex<T>y){return y*=x;}
template<typename T>complex<T>conj(complex<T>x){return x.conj();}
template<typename T>T norm(complex<T>x){return x.norm();}
template<typename T>T abs(complex<T>x){return x.abs();}
template<typename T>T&real(complex<T>&x){return x.real();}
template<typename T>T&imag(complex<T>&x){return x.imag();}
template<typename T>T const real(complex<T>const&x){return x.real();}
template<typename T>T const imag(complex<T>const&x){return x.imag();}
template<typename T>
constexpr complex<T>polar(T r,T theta){return complex<T>::polar(r,theta);}
template<typename T>
std::ostream&operator<<(std::ostream&out,complex<T>x){return out<<x.real()<<' '<<x.imag();}}
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
#line 7 "cp-algo/math/cvector.hpp"
#include <ranges>
#include <bit>
namespace stdx=std::experimental;
namespace cp_algo::math::fft{static constexpr size_t flen=4;
using ftype=double;
using vftype=dx4;
using point=complex<ftype>;
using vpoint=complex<vftype>;
static constexpr vftype vz={};
vpoint vi(vpoint const&r){return{-imag(r),real(r)};}
struct cvector{std::vector<vpoint,big_alloc<vpoint>>r;
cvector(size_t n){n=std::max(flen,std::bit_ceil(n));
r.resize(n/flen);
checkpoint("cvector create");}
vpoint&at(size_t k){return r[k/flen];}
vpoint at(size_t k)const{return r[k/flen];}
template<class pt=point>
void set(size_t k,pt t){if constexpr(std::is_same_v<pt,point>){real(r[k/flen])[k%flen]=real(t);
imag(r[k/flen])[k%flen]=imag(t);}else{at(k)=t;}}
template<class pt=point>
pt get(size_t k)const{if constexpr(std::is_same_v<pt,point>){return{real(r[k/flen])[k%flen],imag(r[k/flen])[k%flen]};}else{return at(k);}}
size_t size()const{return flen*r.size();}
static constexpr size_t eval_arg(size_t n){if(n<pre_evals){return eval_args[n];}else{return eval_arg(n/2)|(n&1)<<(std::bit_width(n)-1);}}
static constexpr point eval_point(size_t n){if(n%2){return-eval_point(n-1);}else if(n%4){return eval_point(n-2)*point(0,1);}else if(n/4<pre_evals){return evalp[n/4];}else{return polar<ftype>(1.,std::numbers::pi/(ftype)std::bit_floor(n)*(ftype)eval_arg(n));}}
static constexpr std::array<point,32>roots=[](){std::array<point,32>res;
for(size_t i=2;i<32;i++){res[i]=polar<ftype>(1.,std::numbers::pi/(1ull<<(i-2)));}
return res;}();
static constexpr point root(size_t n){return roots[std::bit_width(n)];}
template<int step>
static void exec_on_eval(size_t n,size_t k,auto&&callback){callback(k,root(4*step*n)*eval_point(step*k));}
template<int step>
static void exec_on_evals(size_t n,auto&&callback){point factor=root(4*step*n);
for(size_t i=0;i<n;i++){callback(i,factor*eval_point(step*i));}}
void dot(cvector const&t){size_t n=this->size();
exec_on_evals<1>(n/flen,[&](size_t k,point rt){k*=flen;
auto[Ax,Ay]=at(k);
auto Bv=t.at(k);
vpoint res=vz;
for(size_t i=0;i<flen;i++){res+=vpoint(vz+Ax[i],vz+Ay[i])*Bv;
real(Bv)=rotate_right(real(Bv));
imag(Bv)=rotate_right(imag(Bv));
auto x=real(Bv)[0],y=imag(Bv)[0];
real(Bv)[0]=x*real(rt)-y*imag(rt);
imag(Bv)[0]=x*imag(rt)+y*real(rt);}
set(k,res);});
checkpoint("dot");}
template<bool partial=true>
void ifft(){size_t n=size();
if constexpr(!partial){point pi(0,1);
exec_on_evals<4>(n/4,[&](size_t k,point rt){k*=4;
point v1=conj(rt);
point v2=v1*v1;
point v3=v1*v2;
auto A=get(k);
auto B=get(k+1);
auto C=get(k+2);
auto D=get(k+3);
set(k,(A+B)+(C+D));
set(k+2,((A+B)-(C+D))*v2);
set(k+1,((A-B)-pi*(C-D))*v1);
set(k+3,((A-B)+pi*(C-D))*v3);});}
bool parity=std::countr_zero(n)%2;
if(parity){exec_on_evals<2>(n/(2*flen),[&](size_t k,point rt){k*=2*flen;
vpoint cvrt={vz+real(rt),vz-imag(rt)};
auto B=at(k)-at(k+flen);
at(k)+=at(k+flen);
at(k+flen)=B*cvrt;});}
for(size_t leaf=3*flen;leaf<n;leaf+=4*flen){size_t level=std::countr_one(leaf+3);
for(size_t lvl=4+parity;lvl<=level;lvl+=2){size_t i=(1<<lvl)/4;
exec_on_eval<4>(n>>lvl,leaf>>lvl,[&](size_t k,point rt){k<<=lvl;
vpoint v1={vz+real(rt),vz-imag(rt)};
vpoint v2=v1*v1;
vpoint v3=v1*v2;
for(size_t j=k;j<k+i;j+=flen){auto A=at(j);
auto B=at(j+i);
auto C=at(j+2*i);
auto D=at(j+3*i);
at(j)=((A+B)+(C+D));
at(j+2*i)=((A+B)-(C+D))*v2;
at(j+i)=((A-B)-vi(C-D))*v1;
at(j+3*i)=((A-B)+vi(C-D))*v3;}});}}
checkpoint("ifft");
for(size_t k=0;k<n;k+=flen){if constexpr(partial){set(k,get<vpoint>(k)/=vz+ftype(n/flen));}else{set(k,get<vpoint>(k)/=vz+ftype(n));}}}
template<bool partial=true>
void fft(){size_t n=size();
bool parity=std::countr_zero(n)%2;
for(size_t leaf=0;leaf<n;leaf+=4*flen){size_t level=std::countr_zero(n+leaf);
level-=level%2!=parity;
for(size_t lvl=level;lvl>=4;lvl-=2){size_t i=(1<<lvl)/4;
exec_on_eval<4>(n>>lvl,leaf>>lvl,[&](size_t k,point rt){k<<=lvl;
vpoint v1={vz+real(rt),vz+imag(rt)};
vpoint v2=v1*v1;
vpoint v3=v1*v2;
for(size_t j=k;j<k+i;j+=flen){auto A=at(j);
auto B=at(j+i)*v1;
auto C=at(j+2*i)*v2;
auto D=at(j+3*i)*v3;
at(j)=(A+C)+(B+D);
at(j+i)=(A+C)-(B+D);
at(j+2*i)=(A-C)+vi(B-D);
at(j+3*i)=(A-C)-vi(B-D);}});}}
if(parity){exec_on_evals<2>(n/(2*flen),[&](size_t k,point rt){k*=2*flen;
vpoint vrt={vz+real(rt),vz+imag(rt)};
auto t=at(k+flen)*vrt;
at(k+flen)=at(k)-t;
at(k)+=t;});}
if constexpr(!partial){point pi(0,1);
exec_on_evals<4>(n/4,[&](size_t k,point rt){k*=4;
point v1=rt;
point v2=v1*v1;
point v3=v1*v2;
auto A=get(k);
auto B=get(k+1)*v1;
auto C=get(k+2)*v2;
auto D=get(k+3)*v3;
set(k,(A+C)+(B+D));
set(k+1,(A+C)-(B+D));
set(k+2,(A-C)+pi*(B-D));
set(k+3,(A-C)-pi*(B-D));});}
checkpoint("fft");}
static constexpr size_t pre_evals=1<<16;
static const std::array<size_t,pre_evals>eval_args;
static const std::array<point,pre_evals>evalp;};
const std::array<size_t,cvector::pre_evals>cvector::eval_args=[](){std::array<size_t,pre_evals>res={};
for(size_t i=1;i<pre_evals;i++){res[i]=res[i>>1]|(i&1)<<(std::bit_width(i)-1);}
return res;}();
const std::array<point,cvector::pre_evals>cvector::evalp=[](){std::array<point,pre_evals>res={};
res[0]=1;
for(size_t n=1;n<pre_evals;n++){res[n]=polar<ftype>(1.,std::numbers::pi*ftype(eval_args[n])/ftype(4*std::bit_floor(n)));}
return res;}();}
#line 9 "cp-algo/math/fft.hpp"
namespace cp_algo::math::fft{template<modint_type base>
struct dft{cvector A,B;
static base factor,ifactor;
using Int2=base::Int2;
static bool _init;
static int split(){static const int splt=int(std::sqrt(base::mod()))+1;
return splt;}
static uint32_t mod,imod;
static void init(){if(!_init){factor=1+random::rng()%(base::mod()-1);
ifactor=base(1)/factor;
mod=base::mod();
imod=-inv2<uint32_t>(base::mod());
_init=true;}}
dft(size_t n):A(n),B(n){init();}
dft(auto const&a,size_t n,bool partial=true):A(n),B(n){init();
base b2x32=bpow(base(2),32);
u64x4 cur={(bpow(factor,1)*b2x32).getr(),
(bpow(factor,2)*b2x32).getr(),
(bpow(factor,3)*b2x32).getr(),
(bpow(factor,4)*b2x32).getr()};
u64x4 step4=u64x4{}+(bpow(factor,4)*b2x32).getr();
u64x4 stepn=u64x4{}+(bpow(factor,n)*b2x32).getr();
for(size_t i=0;i<std::min(n,std::size(a));i+=flen){auto splt=[&](size_t i,auto mul){if(i>=std::size(a)){return std::pair{vftype(),vftype()};}
u64x4 au={i<std::size(a)?a[i].getr():0,
i+1<std::size(a)?a[i+1].getr():0,
i+2<std::size(a)?a[i+2].getr():0,
i+3<std::size(a)?a[i+3].getr():0};
au=montgomery_mul(au,mul,mod,imod);
au=au>=base::mod()?au-base::mod():au;
auto ai=to_double(i64x4(au>=base::mod()/2?au-base::mod():au));
auto quo=round(ai/split());
return std::pair{ai-quo*split(),quo};};
auto[rai,qai]=splt(i,cur);
auto[rani,qani]=splt(n+i,montgomery_mul(cur,stepn,mod,imod));
A.at(i)=vpoint(rai,rani);
B.at(i)=vpoint(qai,qani);
cur=montgomery_mul(cur,step4,mod,imod);}
checkpoint("dft init");
if(n){if(partial){A.fft();
B.fft();}else{A.template fft<false>();
B.template fft<false>();}}}
template<bool overwrite=true,bool partial=true>
void dot(auto const&C,auto const&D,auto&Aout,auto&Bout,auto&Cout)const{cvector::exec_on_evals<1>(A.size()/flen,[&](size_t k,point rt){k*=flen;
vpoint AC,AD,BC,BD;
AC=AD=BC=BD=vz;
auto Cv=C.at(k),Dv=D.at(k);
if constexpr(partial){auto[Ax,Ay]=A.at(k);
auto[Bx,By]=B.at(k);
for(size_t i=0;i<flen;i++){vpoint Av={vz+Ax[i],vz+Ay[i]},Bv={vz+Bx[i],vz+By[i]};
AC+=Av*Cv;AD+=Av*Dv;
BC+=Bv*Cv;BD+=Bv*Dv;
real(Cv)=rotate_right(real(Cv));
imag(Cv)=rotate_right(imag(Cv));
real(Dv)=rotate_right(real(Dv));
imag(Dv)=rotate_right(imag(Dv));
auto cx=real(Cv)[0],cy=imag(Cv)[0];
auto dx=real(Dv)[0],dy=imag(Dv)[0];
real(Cv)[0]=cx*real(rt)-cy*imag(rt);
imag(Cv)[0]=cx*imag(rt)+cy*real(rt);
real(Dv)[0]=dx*real(rt)-dy*imag(rt);
imag(Dv)[0]=dx*imag(rt)+dy*real(rt);}}else{AC=A.at(k)*Cv;
AD=A.at(k)*Dv;
BC=B.at(k)*Cv;
BD=B.at(k)*Dv;}
if constexpr(overwrite){Aout.at(k)=AC;
Cout.at(k)=AD+BC;
Bout.at(k)=BD;}else{Aout.at(k)+=AC;
Cout.at(k)+=AD+BC;
Bout.at(k)+=BD;}});
checkpoint("dot");}
void dot(auto&&C,auto const&D){dot(C,D,A,B,C);}
void recover_mod(auto&&C,auto&res,size_t k){size_t check=(k+flen-1)/flen*flen;
assert(res.size()>=check);
size_t n=A.size();
auto const splitsplit=base(split()*split()).getr();
base b2x32=bpow(base(2),32);
base b2x64=bpow(base(2),64);
u64x4 cur={(bpow(ifactor,2)*b2x64).getr(),
(bpow(ifactor,3)*b2x64).getr(),
(bpow(ifactor,4)*b2x64).getr(),
(bpow(ifactor,5)*b2x64).getr()};
u64x4 step4=u64x4{}+(bpow(ifactor,4)*b2x32).getr();
u64x4 stepn=u64x4{}+(bpow(ifactor,n)*b2x32).getr();
for(size_t i=0;i<std::min(n,k);i+=flen){auto[Ax,Ay]=A.at(i);
auto[Bx,By]=B.at(i);
auto[Cx,Cy]=C.at(i);
auto set_i=[&](size_t i,auto A,auto B,auto C,auto mul){auto A0=lround(A),A1=lround(C),A2=lround(B);
auto Ai=A0+A1*split()+A2*splitsplit+uint64_t(base::modmod());
auto Au=montgomery_reduce(u64x4(Ai),mod,imod);
Au=montgomery_mul(Au,mul,mod,imod);
Au=Au>=base::mod()?Au-base::mod():Au;
for(size_t j=0;j<flen;j++){res[i+j].setr(typename base::UInt(Au[j]));}};
set_i(i,Ax,Bx,Cx,cur);
if(i+n<k){set_i(i+n,Ay,By,Cy,montgomery_mul(cur,stepn,mod,imod));}
cur=montgomery_mul(cur,step4,mod,imod);}
checkpoint("recover mod");}
void mul(auto&&C,auto const&D,auto&res,size_t k){assert(A.size()==C.size());
size_t n=A.size();
if(!n){res={};
return;}
dot(C,D);
A.ifft();
B.ifft();
C.ifft();
recover_mod(C,res,k);}
void mul_inplace(auto&&B,auto&res,size_t k){mul(B.A,B.B,res,k);}
void mul(auto const&B,auto&res,size_t k){mul(cvector(B.A),B.B,res,k);}
std::vector<base,big_alloc<base>>operator*=(dft&B){std::vector<base,big_alloc<base>>res(2*A.size());
mul_inplace(B,res,2*A.size());
return res;}
std::vector<base,big_alloc<base>>operator*=(dft const&B){std::vector<base,big_alloc<base>>res(2*A.size());
mul(B,res,2*A.size());
return res;}
auto operator*(dft const&B)const{return dft(*this)*=B;}
point operator[](int i)const{return A.get(i);}};
template<modint_type base>base dft<base>::factor=1;
template<modint_type base>base dft<base>::ifactor=1;
template<modint_type base>bool dft<base>::_init=false;
template<modint_type base>uint32_t dft<base>::mod={};
template<modint_type base>uint32_t dft<base>::imod={};
void mul_slow(auto&a,auto const&b,size_t k){if(std::empty(a)||std::empty(b)){a.clear();}else{size_t n=std::min(k,std::size(a));
size_t m=std::min(k,std::size(b));
a.resize(k);
for(int j=int(k-1);j>=0;j--){a[j]*=b[0];
for(int i=std::max(j-(int)n,0)+1;i<std::min(j+1,(int)m);i++){a[j]+=a[j-i]*b[i];}}}}
size_t com_size(size_t as,size_t bs){if(!as||!bs){return 0;}
return std::max(flen,std::bit_ceil(as+bs-1)/2);}
void mul_truncate(auto&a,auto const&b,size_t k){using base=std::decay_t<decltype(a[0])>;
if(std::min({k,std::size(a),std::size(b)})<magic){mul_slow(a,b,k);
return;}
auto n=std::max(flen,std::bit_ceil(
std::min(k,std::size(a))+std::min(k,std::size(b))-1
)/2);
auto A=dft<base>(a|std::views::take(k),n);
auto B=dft<base>(b|std::views::take(k),n);
a.resize((k+flen-1)/flen*flen);
A.mul_inplace(B,a,k);
a.resize(k);}
void mod_split(auto&&x,size_t n,auto k){using base=std::decay_t<decltype(k)>;
dft<base>::init();
assert(std::size(x)==2*n);
u64x4 cur=u64x4{}+(k*bpow(base(2),32)).getr();
for(size_t i=0;i<n;i+=flen){u64x4 xl={x[i].getr(),
x[i+1].getr(),
x[i+2].getr(),
x[i+3].getr()};
u64x4 xr={x[n+i].getr(),
x[n+i+1].getr(),
x[n+i+2].getr(),
x[n+i+3].getr()};
xr=montgomery_mul(xr,cur,dft<base>::mod,dft<base>::imod);
xr=xr>=base::mod()?xr-base::mod():xr;
auto t=xr;
xr=xl-t;
xl+=t;
xl=xl>=base::mod()?xl-base::mod():xl;
xr=xr>=base::mod()?xr+base::mod():xr;
for(size_t k=0;k<flen;k++){x[i+k].setr(typename base::UInt(xl[k]));
x[n+i+k].setr(typename base::UInt(xr[k]));}}
cp_algo::checkpoint("mod split");}
void cyclic_mul(auto&a,auto&&b,size_t k){assert(std::popcount(k)==1);
assert(std::size(a)==std::size(b)&&std::size(a)==k);
using base=std::decay_t<decltype(a[0])>;
dft<base>::init();
if(k<=(1<<16)){std::vector<base,big_alloc<base>>ap(begin(a),end(a));
mul_truncate(ap,b,2*k);
mod_split(ap,k,bpow(dft<base>::factor,k));
std::ranges::copy(ap|std::views::take(k),begin(a));
return;}
k/=2;
auto factor=bpow(dft<base>::factor,k);
mod_split(a,k,factor);
mod_split(b,k,factor);
auto la=std::span(a).first(k);
auto lb=std::span(b).first(k);
auto ra=std::span(a).last(k);
auto rb=std::span(b).last(k);
cyclic_mul(la,lb,k);
auto A=dft<base>(ra,k/2);
auto B=dft<base>(rb,k/2);
A.mul_inplace(B,ra,k);
base i2=base(2).inv();
factor=factor.inv()*i2;
for(size_t i=0;i<k;i++){auto t=(a[i]+a[i+k])*i2;
a[i+k]=(a[i]-a[i+k])*factor;
a[i]=t;}
cp_algo::checkpoint("mod join");}
auto make_copy(auto&&x){return x;}
void cyclic_mul(auto&a,auto const&b,size_t k){return cyclic_mul(a,make_copy(b),k);}
void mul(auto&a,auto&&b){size_t N=size(a)+size(b);
if(N>(1<<20)){N--;
size_t NN=std::bit_ceil(N);
a.resize(NN);
b.resize(NN);
cyclic_mul(a,b,NN);
a.resize(N);}else{mul_truncate(a,b,N-1);}}
void mul(auto&a,auto const&b){size_t N=size(a)+size(b);
if(N>(1<<20)){mul(a,make_copy(b));}else{mul_truncate(a,b,N-1);}}}
#line 12 "cp-algo/math/poly/impl/euclid.hpp"
namespace cp_algo::math::poly::impl{template<typename poly>
using gcd_result=std::pair<
std::list<std::decay_t<poly>>,
linfrac<std::decay_t<poly>>>;
template<typename poly>
gcd_result<poly>half_gcd(poly&&A,poly&&B){assert(A.deg()>=B.deg());
size_t m=size(A.a)/2;
if(B.deg()<(int)m){return{};}
auto[ai,R]=A.divmod(B);
std::tie(A,B)={B,R};
std::list a={ai};
auto T=-linfrac(ai).adj();
auto advance=[&](size_t k){auto[ak,Tk]=half_gcd(A.div_xk(k),B.div_xk(k));
a.splice(end(a),ak);
T.prepend(Tk);
return Tk;};
advance(m).apply(A,B);
if constexpr(std::is_reference_v<poly>){advance(2*m-A.deg()).apply(A,B);}else{advance(2*m-A.deg());}
return{std::move(a),std::move(T)};}
template<typename poly>
gcd_result<poly>full_gcd(poly&&A,poly&&B){using poly_t=std::decay_t<poly>;
std::list<poly_t>ak;
std::vector<linfrac<poly_t>>trs;
while(!B.is_zero()){auto[a0,R]=A.divmod(B);
ak.push_back(a0);
trs.push_back(-linfrac(a0).adj());
std::tie(A,B)={B,R};
auto[a,Tr]=half_gcd(A,B);
ak.splice(end(ak),a);
trs.push_back(Tr);}
return{ak,std::accumulate(rbegin(trs),rend(trs),linfrac<poly_t>{},std::multiplies{})};}
auto convergent(auto L,auto R){using poly=decltype(L)::value_type;
if(R==next(L)){return linfrac(*L);}else{int s=std::transform_reduce(L,R,0,std::plus{},std::mem_fn(&poly::deg));
auto M=L;
for(int c=M->deg();2*c<=s;M++){c+=next(M)->deg();}
return convergent(L,M)*convergent(M,R);}}
template<typename poly>
poly min_rec(poly const&p,size_t d){auto R2=p.mod_xk(d).reversed(d),R1=poly::xk(d);
if(R2.is_zero()){return poly(1);}
auto[a,Tr]=full_gcd(R1,R2);
a.emplace_back();
auto pref=begin(a);
for(int delta=(int)d-a.front().deg();delta>=0;pref++){delta-=pref->deg()+next(pref)->deg();}
return convergent(begin(a),pref).a;}
template<typename poly>
std::optional<poly>inv_mod(poly p,poly q){assert(!q.is_zero());
auto[a,Tr]=full_gcd(q,p);
if(q.deg()!=0){return std::nullopt;}
return Tr.b/q[0];}}
#line 1 "cp-algo/math/poly/impl/div.hpp"
#line 6 "cp-algo/math/poly/impl/div.hpp"
namespace cp_algo::math::poly::impl{auto divmod_slow(auto const&p,auto const&q){auto R=p;
auto D=decltype(p){};
auto q_lead_inv=q.lead().inv();
while(R.deg()>=q.deg()){D.a.push_back(R.lead()*q_lead_inv);
if(D.lead()!=0){for(size_t i=1;i<=q.a.size();i++){R.a[R.a.size()-i]-=D.lead()*q.a[q.a.size()-i];}}
R.a.pop_back();}
std::ranges::reverse(D.a);
R.normalize();
return std::array{D,R};}
template<typename poly>
auto divmod_hint(poly const&p,poly const&q,poly const&qri){assert(!q.is_zero());
int d=p.deg()-q.deg();
if(std::min(d,q.deg())<magic){return divmod_slow(p,q);}
poly D;
if(d>=0){D=(p.reversed().mod_xk(d+1)*qri.mod_xk(d+1)).mod_xk(d+1).reversed(d+1);}
return std::array{D,p-D*q};}
auto divmod(auto const&p,auto const&q){assert(!q.is_zero());
int d=p.deg()-q.deg();
if(std::min(d,q.deg())<magic){return divmod_slow(p,q);}
return divmod_hint(p,q,q.reversed().inv(d+1));}
template<typename poly>
poly powmod_hint(poly const&p,int64_t k,poly const&md,poly const&mdri){return bpow(p,k,poly(1),[&](auto const&p,auto const&q){return divmod_hint(p*q,md,mdri)[1];});}
template<typename poly>
auto powmod(poly const&p,int64_t k,poly const&md){int d=md.deg();
if(p==poly::xk(1)&&false){if(k<md.deg()){return poly::xk(k);}else{auto mdr=md.reversed();
return(mdr.inv(k-md.deg()+1,md.deg())*mdr).reversed(md.deg());}}
if(md==poly::xk(d)){return p.pow(k,d);}
if(md==poly::xk(d)-poly(1)){return p.powmod_circular(k,d);}
return powmod_hint(p,k,md,md.reversed().inv(md.deg()+1));}
template<typename poly>
poly&inv_inplace(poly&q,int64_t k,size_t n){using poly_t=std::decay_t<poly>;
using base=poly_t::base;
if(k<=std::max<int64_t>(n,size(q.a))){return q.inv_inplace(k+n).div_xk_inplace(k);}
if(k%2){return inv_inplace(q,k-1,n+1).div_xk_inplace(1);}
auto[q0,q1]=q.bisect();
auto qq=q0*q0-(q1*q1).mul_xk_inplace(1);
inv_inplace(qq,k/2-q.deg()/2,(n+1)/2+q.deg()/2);
size_t N=fft::com_size(size(q0.a),size(qq.a));
auto q0f=fft::dft<base>(q0.a,N);
auto q1f=fft::dft<base>(q1.a,N);
auto qqf=fft::dft<base>(qq.a,N);
size_t M=q0.deg()+(n+1)/2;
typename poly::Vector A,B;
A.resize((M+fft::flen-1)/fft::flen*fft::flen);
B.resize((M+fft::flen-1)/fft::flen*fft::flen);
q0f.mul(qqf,A,M);
q1f.mul_inplace(qqf,B,M);
q.a.resize(n+1);
for(size_t i=0;i<n;i+=2){q.a[i]=A[q0.deg()+i/2];
q.a[i+1]=-B[q0.deg()+i/2];}
q.a.pop_back();
q.normalize();
return q;}
template<typename poly>
poly&inv_inplace(poly&p,size_t n){using poly_t=std::decay_t<poly>;
using base=poly_t::base;
if(n==1){return p=base(1)/p[0];}
auto[q0,q1]=p.bisect(n);
size_t N=fft::com_size((n+1)/2,(n+1)/2);
auto q0f=fft::dft<base>(q0.a,N);
auto q1f=fft::dft<base>(q1.a,N);
auto qq=poly_t(q0f*q0f)-poly_t(q1f*q1f).mul_xk_inplace(1);
inv_inplace(qq,(n+1)/2);
auto qqf=fft::dft<base>(qq.a,N);
typename poly::Vector A,B;
A.resize(((n+1)/2+fft::flen-1)/fft::flen*fft::flen);
B.resize(((n+1)/2+fft::flen-1)/fft::flen*fft::flen);
q0f.mul(qqf,A,(n+1)/2);
q1f.mul_inplace(qqf,B,(n+1)/2);
p.a.resize(n+1);
for(size_t i=0;i<n;i+=2){p.a[i]=A[i/2];
p.a[i+1]=-B[i/2];}
p.a.pop_back();
p.normalize();
return p;}}
#line 1 "cp-algo/math/combinatorics.hpp"
#line 6 "cp-algo/math/combinatorics.hpp"
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
#line 1 "cp-algo/number_theory/discrete_sqrt.hpp"
#line 6 "cp-algo/number_theory/discrete_sqrt.hpp"
namespace cp_algo::math{template<modint_type base>
std::optional<base>sqrt(base b){if(b==base(0)){return base(0);}else if(bpow(b,(b.mod()-1)/2)!=base(1)){return std::nullopt;}else{while(true){base z=random::rng();
if(z*z==b){return z;}
lin<base>x(1,z,b);
x=bpow(x,(b.mod()-1)/2,lin<base>(0,1,b));
if(x.a!=base(0)){return x.a.inv();}}}}}
#line 15 "cp-algo/math/poly.hpp"
namespace cp_algo::math{template<typename T,class Alloc=big_alloc<T>>
struct poly_t{using Vector=std::vector<T,Alloc>;
using base=T;
Vector a;
poly_t&normalize(){while(deg()>=0&&lead()==base(0)){a.pop_back();}
return*this;}
poly_t(){}
poly_t(T a0):a{a0}{normalize();}
poly_t(Vector const&t):a(t){normalize();}
poly_t(Vector&&t):a(std::move(t)){normalize();}
poly_t&negate_inplace(){std::ranges::transform(a,begin(a),std::negate{});
return*this;}
poly_t operator-()const{return poly_t(*this).negate_inplace();}
poly_t&operator+=(poly_t const&t){a.resize(std::max(size(a),size(t.a)));
std::ranges::transform(a,t.a,begin(a),std::plus{});
return normalize();}
poly_t&operator-=(poly_t const&t){a.resize(std::max(size(a),size(t.a)));
std::ranges::transform(a,t.a,begin(a),std::minus{});
return normalize();}
poly_t operator+(poly_t const&t)const{return poly_t(*this)+=t;}
poly_t operator-(poly_t const&t)const{return poly_t(*this)-=t;}
poly_t&mod_xk_inplace(size_t k){a.resize(std::min(size(a),k));
return normalize();}
poly_t&mul_xk_inplace(size_t k){a.insert(begin(a),k,T(0));
return normalize();}
poly_t&div_xk_inplace(int64_t k){if(k<0){return mul_xk_inplace(-k);}
a.erase(begin(a),begin(a)+std::min<size_t>(k,size(a)));
return normalize();}
poly_t&substr_inplace(size_t l,size_t k){return mod_xk_inplace(l+k).div_xk_inplace(l);}
poly_t mod_xk(size_t k)const{return poly_t(*this).mod_xk_inplace(k);}
poly_t mul_xk(size_t k)const{return poly_t(*this).mul_xk_inplace(k);}
poly_t div_xk(int64_t k)const{return poly_t(*this).div_xk_inplace(k);}
poly_t substr(size_t l,size_t k)const{return poly_t(*this).substr_inplace(l,k);}
poly_t&operator*=(const poly_t&t){fft::mul(a,t.a);normalize();return*this;}
poly_t operator*(const poly_t&t)const{return poly_t(*this)*=t;}
poly_t&operator/=(const poly_t&t){return*this=divmod(t)[0];}
poly_t&operator%=(const poly_t&t){return*this=divmod(t)[1];}
poly_t operator/(poly_t const&t)const{return poly_t(*this)/=t;}
poly_t operator%(poly_t const&t)const{return poly_t(*this)%=t;}
poly_t&operator*=(T const&x){for(auto&it:a){it*=x;}
return normalize();}
poly_t&operator/=(T const&x){return*this*=x.inv();}
poly_t operator*(T const&x)const{return poly_t(*this)*=x;}
poly_t operator/(T const&x)const{return poly_t(*this)/=x;}
poly_t&reverse(size_t n){a.resize(n);
std::ranges::reverse(a);
return normalize();}
poly_t&reverse(){return reverse(size(a));}
poly_t reversed(size_t n)const{return poly_t(*this).reverse(n);}
poly_t reversed()const{return poly_t(*this).reverse();}
std::array<poly_t,2>divmod(poly_t const&b)const{return poly::impl::divmod(*this,b);}
static std::pair<std::list<poly_t>,linfrac<poly_t>>half_gcd(auto&&A,auto&&B){return poly::impl::half_gcd(A,B);}
static std::pair<std::list<poly_t>,linfrac<poly_t>>full_gcd(auto&&A,auto&&B){return poly::impl::full_gcd(A,B);}
static poly_t gcd(poly_t&&A,poly_t&&B){full_gcd(A,B);
return A;}
poly_t min_rec(size_t d)const{return poly::impl::min_rec(*this,d);}
std::optional<poly_t>inv_mod(poly_t const&t)const{return poly::impl::inv_mod(*this,t);};
poly_t negx()const{auto res=*this;
for(int i=1;i<=deg();i+=2){res.a[i]=-res[i];}
return res;}
void print(int n)const{for(int i=0;i<n;i++){std::cout<<(*this)[i]<<' ';}
std::cout<<"\n";}
void print()const{print(deg()+1);}
T eval(T x)const{T res(0);
for(int i=deg();i>=0;i--){res*=x;
res+=a[i];}
return res;}
T lead()const{assert(!is_zero());
return a.back();}
int deg()const{return(int)a.size()-1;}
bool is_zero()const{return a.empty();}
T operator[](int idx)const{return idx<0||idx>deg()?T(0):a[idx];}
T&coef(size_t idx){return a[idx];}
bool operator==(const poly_t&t)const{return a==t.a;}
bool operator!=(const poly_t&t)const{return a!=t.a;}
poly_t&deriv_inplace(int k=1){if(deg()+1<k){return*this=poly_t{};}
for(int i=k;i<=deg();i++){a[i-k]=fact<T>(i)*rfact<T>(i-k)*a[i];}
a.resize(deg()+1-k);
return*this;}
poly_t deriv(int k=1)const{return poly_t(*this).deriv_inplace(k);}
poly_t&integr_inplace(){a.push_back(0);
for(int i=deg()-1;i>=0;i--){a[i+1]=a[i]*small_inv<T>(i+1);}
a[0]=0;
return*this;}
poly_t integr()const{Vector res(deg()+2);
for(int i=0;i<=deg();i++){res[i+1]=a[i]*small_inv<T>(i+1);}
return res;}
size_t trailing_xk()const{if(is_zero()){return-1;}
int res=0;
while(a[res]==T(0)){res++;}
return res;}
poly_t&log_inplace(size_t n){assert(a[0]==T(1));
mod_xk_inplace(n);
return(inv_inplace(n)*=mod_xk_inplace(n).deriv()).mod_xk_inplace(n-1).integr_inplace();}
poly_t log(size_t n)const{return poly_t(*this).log_inplace(n);}
poly_t&mul_truncate(poly_t const&t,size_t k){fft::mul_truncate(a,t.a,k);
return normalize();}
poly_t&exp_inplace(size_t n){if(is_zero()){return*this=T(1);}
assert(a[0]==T(0));
a[0]=1;
size_t a=1;
while(a<n){poly_t C=log(2*a).div_xk_inplace(a)-substr(a,2*a);
*this-=C.mul_truncate(*this,a).mul_xk_inplace(a);
a*=2;}
return mod_xk_inplace(n);}
poly_t exp(size_t n)const{return poly_t(*this).exp_inplace(n);}
poly_t pow_bin(int64_t k,size_t n)const{if(k==0){return poly_t(1).mod_xk(n);}else{auto t=pow(k/2,n);
t=(t*t).mod_xk(n);
return(k%2?*this*t:t).mod_xk(n);}}
poly_t circular_closure(size_t m)const{if(deg()==-1){return*this;}
auto t=*this;
for(size_t i=t.deg();i>=m;i--){t.a[i-m]+=t.a[i];}
t.a.resize(std::min(t.a.size(),m));
return t;}
static poly_t mul_circular(poly_t const&a,poly_t const&b,size_t m){return(a.circular_closure(m)*b.circular_closure(m)).circular_closure(m);}
poly_t powmod_circular(int64_t k,size_t m)const{if(k==0){return poly_t(1);}else{auto t=powmod_circular(k/2,m);
t=mul_circular(t,t,m);
if(k%2){t=mul_circular(t,*this,m);}
return t;}}
poly_t powmod(int64_t k,poly_t const&md)const{return poly::impl::powmod(*this,k,md);}
poly_t pow_dn(int64_t k,size_t n)const{if(n==0){return poly_t(T(0));}
assert((*this)[0]!=T(0));
Vector Q(n);
Q[0]=bpow(a[0],k);
auto a0inv=a[0].inv();
for(int i=1;i<(int)n;i++){for(int j=1;j<=std::min(deg(),i);j++){Q[i]+=a[j]*Q[i-j]*(T(k)*T(j)-T(i-j));}
Q[i]*=small_inv<T>(i)*a0inv;}
return Q;}
poly_t pow(int64_t k,size_t n)const{if(is_zero()){return k?*this:poly_t(1);}
size_t i=trailing_xk();
if(i>0){return k>=int64_t(n+i-1)/(int64_t)i?poly_t(T(0)):div_xk(i).pow(k,n-i*k).mul_xk(i*k);}
if(std::min(deg(),(int)n)<=magic){return pow_dn(k,n);}
if(k<=magic){return pow_bin(k,n);}
T j=a[i];
poly_t t=*this/j;
return bpow(j,k)*(t.log(n)*T(k)).exp(n).mod_xk(n);}
std::optional<poly_t>sqrt(size_t n)const{if(is_zero()){return*this;}
size_t i=trailing_xk();
if(i%2){return std::nullopt;}else if(i>0){auto ans=div_xk(i).sqrt(n-i/2);
return ans?ans->mul_xk(i/2):ans;}
auto st=math::sqrt((*this)[0]);
if(st){poly_t ans=*st;
size_t a=1;
while(a<n){a*=2;
ans-=(ans-mod_xk(a)*ans.inv(a)).mod_xk(a)/2;}
return ans.mod_xk(n);}
return std::nullopt;}
poly_t mulx(T a)const{T cur=1;
poly_t res(*this);
for(int i=0;i<=deg();i++){res.coef(i)*=cur;
cur*=a;}
return res;}
poly_t mulx_sq(T a)const{T cur=1,total=1;
poly_t res(*this);
for(int i=0;i<=deg();i++){res.coef(i)*=total;
cur*=a;
total*=cur;}
return res;}
poly_t chirpz(T z,int n)const{if(is_zero()){return Vector(n,0);}
if(z==T(0)){Vector ans(n,(*this)[0]);
if(n>0){ans[0]=accumulate(begin(a),end(a),T(0));}
return ans;}
auto A=mulx_sq(z.inv());
auto B=ones(n+deg()).mulx_sq(z);
return semicorr(B,A).mod_xk(n).mulx_sq(z.inv());}
static auto _1mzk_prod_inv(T z,int n){Vector res(n,1),zk(n);
zk[0]=1;
for(int i=1;i<n;i++){zk[i]=zk[i-1]*z;
res[i]=res[i-1]*(T(1)-zk[i]);}
res.back()=res.back().inv();
for(int i=n-2;i>=0;i--){res[i]=(T(1)-zk[i+1])*res[i+1];}
return res;}
static auto _1mzkx_prod(T z,int n){if(n==1){return poly_t(Vector{1,-1});}else{auto t=_1mzkx_prod(z,n/2);
t*=t.mulx(bpow(z,n/2));
if(n%2){t*=poly_t(Vector{1,-bpow(z,n-1)});}
return t;}}
poly_t chirpz_inverse(T z,int n)const{if(is_zero()){return{};}
if(z==T(0)){if(n==1){return*this;}else{return Vector{(*this)[1],(*this)[0]-(*this)[1]};}}
Vector y(n);
for(int i=0;i<n;i++){y[i]=(*this)[i];}
auto prods_pos=_1mzk_prod_inv(z,n);
auto prods_neg=_1mzk_prod_inv(z.inv(),n);
T zn=bpow(z,n-1).inv();
T znk=1;
for(int i=0;i<n;i++){y[i]*=znk*prods_neg[i]*prods_pos[(n-1)-i];
znk*=zn;}
poly_t p_over_q=poly_t(y).chirpz(z,n);
poly_t q=_1mzkx_prod(z,n);
return(p_over_q*q).mod_xk_inplace(n).reverse(n);}
static poly_t build(std::vector<poly_t>&res,int v,auto L,auto R){if(R-L==1){return res[v]=Vector{-*L,1};}else{auto M=L+(R-L)/2;
return res[v]=build(res,2*v,L,M)*build(res,2*v+1,M,R);}}
poly_t to_newton(std::vector<poly_t>&tree,int v,auto l,auto r){if(r-l==1){return*this;}else{auto m=l+(r-l)/2;
auto A=(*this%tree[2*v]).to_newton(tree,2*v,l,m);
auto B=(*this/tree[2*v]).to_newton(tree,2*v+1,m,r);
return A+B.mul_xk(m-l);}}
poly_t to_newton(Vector p){if(is_zero()){return*this;}
size_t n=p.size();
std::vector<poly_t>tree(4*n);
build(tree,1,begin(p),end(p));
return to_newton(tree,1,begin(p),end(p));}
Vector eval(std::vector<poly_t>&tree,int v,auto l,auto r){if(r-l==1){return{eval(*l)};}else{auto m=l+(r-l)/2;
auto A=(*this%tree[2*v]).eval(tree,2*v,l,m);
auto B=(*this%tree[2*v+1]).eval(tree,2*v+1,m,r);
A.insert(end(A),begin(B),end(B));
return A;}}
Vector eval(Vector x){size_t n=x.size();
if(is_zero()){return Vector(n,T(0));}
std::vector<poly_t>tree(4*n);
build(tree,1,begin(x),end(x));
return eval(tree,1,begin(x),end(x));}
poly_t inter(std::vector<poly_t>&tree,int v,auto ly,auto ry){if(ry-ly==1){return{*ly/a[0]};}else{auto my=ly+(ry-ly)/2;
auto A=(*this%tree[2*v]).inter(tree,2*v,ly,my);
auto B=(*this%tree[2*v+1]).inter(tree,2*v+1,my,ry);
return A*tree[2*v+1]+B*tree[2*v];}}
static auto inter(Vector x,Vector y){size_t n=x.size();
std::vector<poly_t>tree(4*n);
return build(tree,1,begin(x),end(x)).deriv().inter(tree,1,begin(y),end(y));}
static auto resultant(poly_t a,poly_t b){if(b.is_zero()){return 0;}else if(b.deg()==0){return bpow(b.lead(),a.deg());}else{int pw=a.deg();
a%=b;
pw-=a.deg();
auto mul=bpow(b.lead(),pw)*T((b.deg()&a.deg()&1)?-1:1);
auto ans=resultant(b,a);
return ans*mul;}}
static poly_t xk(size_t n){return poly_t(T(1)).mul_xk(n);}
static poly_t ones(size_t n){return Vector(n,1);}
static poly_t expx(size_t n){return ones(n).borel();}
static poly_t log1px(size_t n){Vector coeffs(n,0);
for(size_t i=1;i<n;i++){coeffs[i]=(i&1?T(i).inv():-T(i).inv());}
return coeffs;}
static poly_t log1mx(size_t n){return-ones(n).integr();}
static poly_t corr(poly_t const&a,poly_t const&b){return a*b.reversed();}
static poly_t semicorr(poly_t const&a,poly_t const&b){return corr(a,b).div_xk(b.deg());}
poly_t invborel()const{auto res=*this;
for(int i=0;i<=deg();i++){res.coef(i)*=fact<T>(i);}
return res;}
poly_t borel()const{auto res=*this;
for(int i=0;i<=deg();i++){res.coef(i)*=rfact<T>(i);}
return res;}
poly_t shift(T a)const{return semicorr(invborel(),expx(deg()+1).mulx(a)).borel();}
poly_t x2(){Vector res(2*a.size());
for(size_t i=0;i<a.size();i++){res[2*i]=a[i];}
return res;}
std::array<poly_t,2>bisect(size_t n)const{n=std::min(n,size(a));
Vector res[2];
for(size_t i=0;i<n;i++){res[i%2].push_back(a[i]);}
return{res[0],res[1]};}
std::array<poly_t,2>bisect()const{return bisect(size(a));}
static T kth_rec_inplace(poly_t&P,poly_t&Q,int64_t k){while(k>Q.deg()){size_t n=Q.a.size();
auto[Q0,Q1]=Q.bisect();
auto[P0,P1]=P.bisect();
size_t N=fft::com_size((n+1)/2,(n+1)/2);
auto Q0f=fft::dft<T>(Q0.a,N);
auto Q1f=fft::dft<T>(Q1.a,N);
auto P0f=fft::dft<T>(P0.a,N);
auto P1f=fft::dft<T>(P1.a,N);
Q=poly_t(Q0f*Q0f)-=poly_t(Q1f*Q1f).mul_xk_inplace(1);
if(k%2){P=poly_t(Q0f*=P1f)-=poly_t(Q1f*=P0f);}else{P=poly_t(Q0f*=P0f)-=poly_t(Q1f*=P1f).mul_xk_inplace(1);}
k/=2;}
return(P*=Q.inv_inplace(Q.deg()+1))[(int)k];}
static T kth_rec(poly_t const&P,poly_t const&Q,int64_t k){return kth_rec_inplace(poly_t(P),poly_t(Q),k);}
poly_t&inv_inplace(size_t n){return poly::impl::inv_inplace(*this,n);}
poly_t inv(size_t n)const{return poly_t(*this).inv_inplace(n);}
poly_t&inv_inplace(int64_t k,size_t n){return poly::impl::inv_inplace(*this,k,n);}
poly_t inv(int64_t k,size_t n)const{return poly_t(*this).inv_inplace(k,n);}
static poly_t compose(poly_t A,poly_t B,int n){int q=std::sqrt(n);
std::vector<poly_t>Bk(q);
auto Bq=B.pow(q,n);
Bk[0]=poly_t(T(1));
for(int i=1;i<q;i++){Bk[i]=(Bk[i-1]*B).mod_xk(n);}
poly_t Bqk(1);
poly_t ans;
for(int i=0;i<=n/q;i++){poly_t cur;
for(int j=0;j<q;j++){cur+=Bk[j]*A[i*q+j];}
ans+=(Bqk*cur).mod_xk(n);
Bqk=(Bqk*Bq).mod_xk(n);}
return ans;}
static poly_t compose_large(poly_t A,poly_t B,int n){if(B[0]!=T(0)){return compose_large(A.shift(B[0]),B-B[0],n);}
int q=std::sqrt(n);
auto[B0,B1]=std::make_pair(B.mod_xk(q),B.div_xk(q));
B0=B0.div_xk(1);
std::vector<poly_t>pw(A.deg()+1);
auto getpow=[&](int k){return pw[k].is_zero()?pw[k]=B0.pow(k,n-k):pw[k];};
std::function<poly_t(poly_t const&,int,int)>compose_dac=[&getpow,&compose_dac](poly_t const&f,int m,int N){if(f.deg()<=0){return f;}
int k=m/2;
auto[f0,f1]=std::make_pair(f.mod_xk(k),f.div_xk(k));
auto[A,B]=std::make_pair(compose_dac(f0,k,N),compose_dac(f1,m-k,N-k));
return(A+(B.mod_xk(N-k)*getpow(k).mod_xk(N-k)).mul_xk(k)).mod_xk(N);};
int r=n/q;
auto Ar=A.deriv(r);
auto AB0=compose_dac(Ar,Ar.deg()+1,n);
auto Bd=B0.mul_xk(1).deriv();
poly_t ans=T(0);
std::vector<poly_t>B1p(r+1);
B1p[0]=poly_t(T(1));
for(int i=1;i<=r;i++){B1p[i]=(B1p[i-1]*B1.mod_xk(n-i*q)).mod_xk(n-i*q);}
while(r>=0){ans+=(AB0.mod_xk(n-r*q)*rfact<T>(r)*B1p[r]).mul_xk(r*q).mod_xk(n);
r--;
if(r>=0){AB0=((AB0*Bd).integr()+A[r]*fact<T>(r)).mod_xk(n);}}
return ans;}};
template<typename base>
static auto operator*(const auto&a,const poly_t<base>&b){return b*a;}};
#line 1 "cp-algo/linalg/matrix.hpp"
#line 1 "cp-algo/linalg/vector.hpp"
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
#line 8 "cp-algo/linalg/frobenius.hpp"
namespace cp_algo::linalg{enum frobenius_mode{blocks,full};
template<frobenius_mode mode=blocks>
auto frobenius_form(auto const&A){using matrix=std::decay_t<decltype(A)>;
using vec_t=matrix::vec_t;
using base=typename matrix::base;
using base=matrix::base;
using polyn=math::poly_t<base>;
assert(A.n()==A.m());
size_t n=A.n();
std::vector<polyn>charps;
std::vector<vec_t>basis,basis_init;
while(size(basis)<n){size_t start=size(basis);
auto generate_block=[&](auto x){while(true){vec_t y=x|vec_t::ei(n+1,size(basis));
for(auto&it:basis){y.reduce_by(it);}
y.normalize();
if(std::ranges::count(y|std::views::take(n),base(0))==int(n)){return polyn(typename polyn::Vector(begin(y)+n,end(y)));}else{basis_init.push_back(x);
basis.push_back(y);
x=A.apply(x);}}};
auto full_rec=generate_block(vec_t::random(n));
if constexpr(mode==full){if(full_rec.mod_xk(start)!=polyn()){auto charp=full_rec.div_xk(start);
auto x=basis_init[start];
auto shift=full_rec/charp;
for(int j=0;j<shift.deg();j++){x.add_scaled(basis_init[j],shift[j]);}
basis.resize(start);
basis_init.resize(start);
full_rec=generate_block(x.normalize());}}
charps.push_back(full_rec.div_xk(start));}
if constexpr(mode==full){for(size_t i=0;i<n;i++){for(size_t j=i+1;j<n;j++){basis[i].reduce_by(basis[j]);}
basis[i].normalize();}
auto T=matrix(basis_init);
auto Tinv=matrix(basis);
std::ignore=Tinv.sort_classify(n);
for(size_t i=0;i<n;i++){Tinv[i]=vec_t(
Tinv[i]|std::views::drop(n)|std::views::take(n)
)*(base(1)/Tinv[i][i]);}
return std::tuple{T,Tinv,charps};}else{return charps;}}
template<typename base>
auto with_frobenius(matrix<base>const&A,auto&&callback){auto[T,Tinv,charps]=frobenius_form<full>(A);
std::vector<matrix<base>>blocks;
for(auto charp:charps){matrix<base>block(charp.deg());
auto xk=callback(charp);
for(size_t i=0;i<block.n();i++){std::ranges::copy(xk.a,begin(block[i]));
xk=xk.mul_xk(1)%charp;}
blocks.push_back(block);}
auto S=matrix<base>::block_diagonal(blocks);
return Tinv*S*T;}
template<typename base>
auto frobenius_pow(matrix<base>const&A,uint64_t k){return with_frobenius(A,[k](auto const&charp){return math::poly_t<base>::xk(1).powmod(k,charp);});}};
#line 8 "verify/linalg/characteristic.test.cpp"
using namespace std;
using namespace cp_algo::math;
using namespace cp_algo::linalg;
const int64_t mod=998244353;
using base=modint<mod>;
using polyn=poly_t<base>;
void solve(){size_t n;
cin>>n;
matrix<base>A(n);
A.read();
auto blocks=frobenius_form(A);
reduce(begin(blocks),end(blocks),polyn(1),multiplies{}).print();}
signed main(){ios::sync_with_stdio(0);
cin.tie(0);
int t=1;
while(t--){solve();}}