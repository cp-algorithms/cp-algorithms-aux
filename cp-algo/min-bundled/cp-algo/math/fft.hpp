#line 1 "cp-algo/math/fft.hpp"
#line 1 "cp-algo/number_theory/modint.hpp"
#line 1 "cp-algo/math/common.hpp"
#include <functional>
#include <cstdint>
#include <cassert>
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
#line 4 "cp-algo/number_theory/modint.hpp"
#include <iostream>
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
#line 4 "cp-algo/util/checkpoint.hpp"
#include <chrono>
#include <string>
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
#line 1 "cp-algo/random/rng.hpp"
#line 4 "cp-algo/random/rng.hpp"
#include <random>
namespace cp_algo::random{std::mt19937_64 gen(
std::chrono::steady_clock::now().time_since_epoch().count()
);
uint64_t rng(){return gen();}}
#line 1 "cp-algo/math/cvector.hpp"
#line 1 "cp-algo/util/simd.hpp"
#include <experimental/simd>
#line 5 "cp-algo/util/simd.hpp"
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
#line 1 "cp-algo/util/complex.hpp"
#line 4 "cp-algo/util/complex.hpp"
#include <cmath>
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
#include <vector>
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