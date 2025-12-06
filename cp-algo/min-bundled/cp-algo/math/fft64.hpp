#line 1 "cp-algo/math/fft64.hpp"
#line 1 "cp-algo/random/rng.hpp"
#include <chrono>
#include <random>
namespace cp_algo::random{std::mt19937_64 gen(
std::chrono::steady_clock::now().time_since_epoch().count()
);
uint64_t rng(){return gen();}}
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
#include <iostream>
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
#line 1 "cp-algo/util/checkpoint.hpp"
#line 5 "cp-algo/util/checkpoint.hpp"
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