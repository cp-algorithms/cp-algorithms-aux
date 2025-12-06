#line 1 "verify/poly/convolution64.test.cpp"
#define PROBLEM "https:#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_CHECKPOINT
#include <bits/stdc++.h>
#line 1 "cp-algo/math/fft64.hpp"
#line 1 "cp-algo/random/rng.hpp"
#line 5 "cp-algo/random/rng.hpp"
namespace cp_algo::random{std::mt19937_64 gen(
std::chrono::steady_clock::now().time_since_epoch().count()
);
uint64_t rng(){return gen();}}
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
#line 6 "cp-algo/math/fft64.hpp"
namespace cp_algo::math::fft{struct dft64{std::vector<cp_algo::math::fft::cvector>cv;
static uint64_t factor,ifactor;
static bool _init;
static void init(){if(_init)return;
_init=true;
factor=random::rng();
if(factor%2==0){factor++;}
ifactor=inv2(factor);}
dft64(auto const&a,size_t n):cv(4,n){init();
uint64_t cur=1,step=bpow(factor,n);
for(size_t i=0;i<std::min(std::size(a),n);i++){auto split=[&](size_t i,uint64_t mul)->std::array<int16_t,4>{uint64_t x=i<std::size(a)?a[i]*mul:0;
std::array<int16_t,4>res;
for(int z=0;z<4;z++){res[z]=int16_t(x);
x=(x>>16)+(res[z]<0);}
return res;};
auto re=split(i,cur);
auto im=split(n+i,cur*step);
for(int z=0;z<4;z++){real(cv[z].at(i))[i%4]=re[z];
imag(cv[z].at(i))[i%4]=im[z];}
cur*=factor;}
checkpoint("dft64 init");
for(auto&x:cv){x.fft();}}
void dot(dft64 const&t){size_t N=cv[0].size();
cvector::exec_on_evals<1>(N/flen,[&](size_t k,point rt){k*=flen;
auto[A0x,A0y]=cv[0].at(k);
auto[A1x,A1y]=cv[1].at(k);
auto[A2x,A2y]=cv[2].at(k);
auto[A3x,A3y]=cv[3].at(k);
std::array B={t.cv[0].at(k),
t.cv[1].at(k),
t.cv[2].at(k),
t.cv[3].at(k)};
std::array<vpoint,4>C={vz,vz,vz,vz};
for(size_t i=0;i<flen;i++){std::array A={vpoint{vz+A0x[i],vz+A0y[i]},
vpoint{vz+A1x[i],vz+A1y[i]},
vpoint{vz+A2x[i],vz+A2y[i]},
vpoint{vz+A3x[i],vz+A3y[i]}};
for(size_t k=0;k<4;k++){for(size_t i=0;i<=k;i++){C[k]+=A[i]*B[k-i];}}
for(size_t k=0;k<4;k++){real(B[k])=rotate_right(real(B[k]));
imag(B[k])=rotate_right(imag(B[k]));
auto bx=real(B[k])[0],by=imag(B[k])[0];
real(B[k])[0]=bx*real(rt)-by*imag(rt);
imag(B[k])[0]=bx*imag(rt)+by*real(rt);}}
cv[0].at(k)=C[0];
cv[1].at(k)=C[1];
cv[2].at(k)=C[2];
cv[3].at(k)=C[3];});
checkpoint("dot");
for(auto&x:cv){x.ifft();}}
void recover_mod(auto&res,size_t k){size_t n=cv[0].size();
uint64_t cur=1,step=bpow(ifactor,n);
for(size_t i=0;i<std::min(k,n);i++){std::array re={real(cv[0].get(i)),real(cv[1].get(i)),real(cv[2].get(i)),real(cv[3].get(i))};
std::array im={imag(cv[0].get(i)),imag(cv[1].get(i)),imag(cv[2].get(i)),imag(cv[3].get(i))};
auto set_i=[&](size_t i,auto&x,auto mul){if(i>=k)return;
res[i]=llround(x[0])+(llround(x[1])<<16)+(llround(x[2])<<32)+(llround(x[3])<<48);
res[i]*=mul;};
set_i(i,re,cur);
set_i(n+i,im,cur*step);
cur*=ifactor;}
cp_algo::checkpoint("recover mod");}};
uint64_t dft64::factor=1,dft64::ifactor=1;
bool dft64::_init=false;
void conv64(auto&a,auto const&b){size_t n=a.size(),m=b.size();
size_t N=std::max(flen,std::bit_ceil(n+m-1)/2);
dft64 A(a,N),B(b,N);
A.dot(B);
a.resize(n+m-1);
A.recover_mod(a,n+m-1);}}
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
#line 8 "verify/poly/convolution64.test.cpp"
using namespace std;
void solve(){int n,m;
cin>>n>>m;
vector<uint64_t,cp_algo::big_alloc<uint64_t>>a(n),b(m);
for(auto&x:a)cin>>x;
for(auto&x:b)cin>>x;
cp_algo::checkpoint("read");
cp_algo::math::fft::conv64(a,b);
for(auto x:a){cout<<uint64_t(x)<<" ";}
cp_algo::checkpoint("write");
cp_algo::checkpoint<1>();}
signed main(){ios::sync_with_stdio(0);
cin.tie(0);
solve();}