#line 1 "verify/structures/segtree/range_affine_range_sum.test.cpp"
#define PROBLEM "https:#line 1 "cp-algo/structures/segtree/metas/affine.hpp"
#line 1 "cp-algo/structures/segtree/metas/base.hpp"
#include <cstddef>
namespace cp_algo::structures::segtree::metas{template<typename derived_meta>
struct base_meta{using meta=derived_meta;
virtual void pull(meta const&,meta const&,size_t,size_t){};
virtual void push(meta*,meta*,size_t,size_t){};};}
#line 1 "cp-algo/math/affine.hpp"
#include <optional>
#include <utility>
#include <cassert>
#include <tuple>
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
#line 5 "cp-algo/structures/segtree/metas/affine.hpp"
namespace cp_algo::structures::segtree::metas{template<typename base>
struct affine_meta:base_meta<affine_meta<base>>{using meta=affine_meta;
using lin=math::lin<base>;
base sum=0;
lin to_push={};
affine_meta(){}
affine_meta(base sum):sum(sum){}
void push(meta*L,meta*R,size_t l,size_t r)override{if(to_push.a!=1||to_push.b!=0){sum=to_push.a*sum+to_push.b*(r-l);
if(r-l>1){L->to_push.prepend(to_push);
R->to_push.prepend(to_push);}
to_push={};}}
void pull(meta const&L,meta const&R,size_t,size_t)override{sum=L.sum+R.sum;}};}
#line 1 "cp-algo/structures/segtree.hpp"
#include <vector>
#include <numeric>
namespace cp_algo::structures{template<typename meta>
struct segtree_t{const size_t N;
std::vector<meta>_meta;
segtree_t(size_t n):N(n),_meta(4*N){}
segtree_t(std::vector<meta>leafs):N(size(leafs)),_meta(4*N){build(leafs);}
void pull(size_t v,size_t l,size_t r){if(r-l>1){_meta[v].pull(_meta[2*v],_meta[2*v+1],l,r);}}
void push(size_t v,size_t l,size_t r){if(r-l>1){_meta[v].push(&_meta[2*v],&_meta[2*v+1],l,r);}else{_meta[v].push(nullptr,nullptr,l,r);}}
void build(auto&a,size_t v,size_t l,size_t r){if(r-l==1){if(l<size(a)){_meta[v]=a[l];}}else{size_t m=std::midpoint(l,r);
build(a,2*v,l,m);
build(a,2*v+1,m,r);
pull(v,l,r);}}
void build(auto&a){build(a,1,0,N);}
void exec_on_segment(size_t a,size_t b,auto func,auto proceed,auto stop,size_t v,size_t l,size_t r){push(v,l,r);
if(r<=a||b<=l||stop(_meta[v])){return;}else if(a<=l&&r<=b&&proceed(_meta[v])){func(_meta[v]);
push(v,l,r);}else{size_t m=std::midpoint(l,r);
exec_on_segment(a,b,func,proceed,stop,2*v,l,m);
exec_on_segment(a,b,func,proceed,stop,2*v+1,m,r);
pull(v,l,r);}}
static constexpr auto default_true=[](auto const&){return true;};
static constexpr auto default_false=[](auto const&){return false;};
void exec_on_segment(size_t a,size_t b,auto func,auto proceed,auto stop){exec_on_segment(a,b,func,proceed,stop,1,0,N);}
void exec_on_segment(size_t a,size_t b,auto func){exec_on_segment(a,b,func,default_true,default_false);}};}
#line 1 "cp-algo/number_theory/modint.hpp"
#line 1 "cp-algo/math/common.hpp"
#include <functional>
#include <cstdint>
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
#line 6 "verify/structures/segtree/range_affine_range_sum.test.cpp"
#include <bits/stdc++.h>
using namespace std;
using namespace cp_algo::structures;
using base=cp_algo::math::modint<998244353>;
using meta=segtree::metas::affine_meta<base>;
void solve(){int n,q;
cin>>n>>q;
vector<meta>a(n);
for(int i=0;i<n;i++){int ai;
cin>>ai;
a[i]={ai};}
segtree_t<meta>me(a);
while(q--){int t;
cin>>t;
if(t==0){int l,r,b,c;
cin>>l>>r>>b>>c;
me.exec_on_segment(l,r,[&](auto&meta){meta.to_push.prepend(meta::lin(b,c));});}else{int l,r;
cin>>l>>r;
base ans=0;
me.exec_on_segment(l,r,[&](auto meta){ans+=meta.sum;});
cout<<ans<<"\n";}}}
signed main(){ios::sync_with_stdio(0);
cin.tie(0);
int t=1;
while(t--){solve();}}