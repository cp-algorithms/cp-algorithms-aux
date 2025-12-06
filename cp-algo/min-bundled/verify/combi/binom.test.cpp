#line 1 "verify/combi/binom.test.cpp"
#define PROBLEM "https:#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_MAXN 1e7
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
#line 1 "cp-algo/math/combinatorics.hpp"
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
#line 7 "verify/combi/binom.test.cpp"
#include <bits/stdc++.h>
using namespace std;
using namespace cp_algo;
using namespace math;
using base=dynamic_modint<>;
void solve(){int n,r;
cin>>n>>r;
cout<<binom<base>(n,r)<<"\n";}
signed main(){ios::sync_with_stdio(0);
cin.tie(0);
int t=1;
cin>>t;
int m;
cin>>m;
base::switch_mod(m);
while(t--){solve();}}