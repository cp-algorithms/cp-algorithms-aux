#line 1 "cp-algo/number_theory/two_squares.hpp"
#line 1 "cp-algo/number_theory/euler.hpp"
#line 1 "cp-algo/number_theory/factorize.hpp"
#line 1 "cp-algo/number_theory/primality.hpp"
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
#line 4 "cp-algo/number_theory/primality.hpp"
#include <algorithm>
#include <bit>
namespace cp_algo::math{template<typename _Int>
bool is_prime(_Int m){using Int=std::make_signed_t<_Int>;
using UInt=std::make_unsigned_t<Int>;
if(m==1||m%2==0){return m==2;}
int s=std::countr_zero(UInt(m-1));
auto d=(m-1)>>s;
using base=dynamic_modint<Int>;
auto test=[&](base x){x=bpow(x,d);
if(std::abs(x.rem())<=1){return true;}
for(int i=1;i<s&&x!=-1;i++){x*=x;}
return x==-1;};
return base::with_mod(m,[&](){#ifdef CP_ALGO_NUMBER_THEORY_PRIMALITY_BASES_HPP
uint16_t base2=7,base3=61;
if(m!=uint32_t(m)){base2=base_table1[uint32_t(m*0xAD625B89)>>18];
base3=base_table2[base2>>13];}
return test(2)&&test(base2)&&test(base3);
#else
return std::ranges::all_of(std::array{2,325,9375,28178,450775,9780504,1795265022},test);
#endif});}}
#line 1 "cp-algo/random/rng.hpp"
#include <chrono>
#include <random>
namespace cp_algo::random{std::mt19937_64 gen(
std::chrono::steady_clock::now().time_since_epoch().count()
);
uint64_t rng(){return gen();}}
#line 5 "cp-algo/number_theory/factorize.hpp"
#include <generator>
namespace cp_algo::math{template<typename _Int>
auto proper_divisor(_Int m){using Int=std::make_signed_t<_Int>;
using base=dynamic_modint<Int>;
return m%2==0?2:base::with_mod(m,[&](){base t=random::rng();
auto f=[&](auto x){return x*x+t;};
base x=0,y=0;
base g=1;
while(g==1){for(int i=0;i<64;i++){x=f(x);
y=f(f(y));
if(x==y)[[unlikely]]{t=random::rng();
x=y=0;}else{base t=g*(x-y);
g=t==0?g:t;}}
g=std::gcd(g.getr(),m);}
return g.getr();});}
template<typename Int>
std::generator<Int>factorize(Int m){if(is_prime(m)){co_yield m;}else if(m>1){auto g=proper_divisor(m);
co_yield std::ranges::elements_of(factorize(g));
co_yield std::ranges::elements_of(factorize(m/g));}}
template<typename Int>
std::generator<Int>divisors_sqrt(Int m){for(Int i=1;i*i<=m;i++){if(m%i==0){co_yield i;
if(i*i!=m){co_yield m/i;}}}}}
#line 4 "cp-algo/number_theory/euler.hpp"
namespace cp_algo::math{auto euler_phi(auto m){auto primes=to<std::vector>(factorize(m));
std::ranges::sort(primes);
auto[from,to]=std::ranges::unique(primes);
primes.erase(from,to);
auto ans=m;
for(auto it:primes){ans-=ans/it;}
return ans;}
template<modint_type base>
auto period(base x){auto ans=euler_phi(base::mod());
base x0=bpow(x,ans);
for(auto t:factorize(ans)){while(ans%t==0&&x0*bpow(x,ans/t)==x0){ans/=t;}}
return ans;}
template<typename _Int>
_Int primitive_root(_Int p){using Int=std::make_signed_t<_Int>;
using base=dynamic_modint<Int>;
return base::with_mod(p,[p](){base t=1;
while(period(t)!=p-1){t=random::rng();}
return t.getr();});}}
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
#line 6 "cp-algo/number_theory/two_squares.hpp"
#include <utility>
#include <vector>
#include <map>
namespace cp_algo::math{template<typename T>
using gaussint=complex<T>;
template<typename _Int>
auto two_squares_prime_any(_Int p){if(p==2){return gaussint<_Int>{1,1};}
assert(p%4==1);
using Int=std::make_signed_t<_Int>;
using base=dynamic_modint<Int>;
return base::with_mod(p,[&](){base g=primitive_root(p);
int64_t i=bpow(g,(p-1)/4).getr();
int64_t q0=1,q1=0;
int64_t r=i,m=p;
do{int64_t d=r/m;
q0=std::exchange(q1,q0+d*q1);
r=std::exchange(m,r%m);}while(q1<p/q1);
return gaussint<_Int>{q0,(base(i)*base(q0)).rem()};});}
template<typename Int>
std::vector<gaussint<Int>>two_squares_all(Int n){if(n==0){return{0};}
auto primes=factorize(n);
std::map<Int,int>cnt;
for(auto p:primes){cnt[p]++;}
std::vector<gaussint<Int>>res={1};
for(auto[p,c]:cnt){std::vector<gaussint<Int>>nres;
if(p%4==3){if(c%2==0){auto mul=bpow(gaussint<Int>(p),c/2);
for(auto p:res){nres.push_back(p*mul);}}}else if(p%4==1){auto base=two_squares_prime_any(p);
for(int i=0;i<=c;i++){auto mul=bpow(base,i)*bpow(conj(base),c-i);
for(auto p:res){nres.push_back(p*mul);}}}else if(p%4==2){auto mul=bpow(gaussint<Int>(1,1),c);
for(auto p:res){nres.push_back(p*mul);}}
res=nres;}
std::vector<gaussint<Int>>nres;
for(auto p:res){while(p.real()<0||p.imag()<0){p*=gaussint<Int>(0,1);}
nres.push_back(p);
if(!p.real()||!p.imag()){nres.emplace_back(p.imag(),p.real());}}
return nres;}}