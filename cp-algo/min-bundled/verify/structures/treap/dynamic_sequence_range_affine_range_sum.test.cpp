#line 1 "verify/structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp"
#define PROBLEM "https:#line 1 "cp-algo/number_theory/modint.hpp"
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
#line 1 "cp-algo/structures/treap/metas/reverse.hpp"
#line 1 "cp-algo/structures/treap/metas/base.hpp"
#line 1 "cp-algo/structures/treap/common.hpp"
#define _safe(t, op) (t ? t->op : typename std::remove_reference_t<decltype(t->op)>())
#line 5 "cp-algo/structures/treap/metas/base.hpp"
#include <algorithm>
#line 7 "cp-algo/structures/treap/metas/base.hpp"
#define _safe_meta(i, op) _safe(i, _meta.op)
namespace cp_algo::structures::treap::metas{struct base_meta{void pull(auto const,auto const){}
void push(auto&,auto&){}};}
#line 1 "cp-algo/math/affine.hpp"
#include <optional>
#include <utility>
#line 6 "cp-algo/math/affine.hpp"
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
#line 6 "cp-algo/structures/treap/metas/reverse.hpp"
namespace cp_algo::structures::treap::metas{template<typename base>
struct reverse_meta:base_meta{using lin=math::lin<base>;
base val;
size_t sz=1;
bool reverse=false;
base sum=val;
lin to_push={};
reverse_meta(base val):val(val){}
void pull(auto const L,auto const R){sum=val+_safe_meta(L,sum)+_safe_meta(R,sum);
sz=1+_safe_meta(L,sz)+_safe_meta(R,sz);}
void add_push(lin const&t){val=t.eval(val);
sum=t.a*sum+t.b*sz;
to_push.prepend(t);}
void push(auto&L,auto&R){if(reverse){reverse=false;
std::swap(L,R);
_safe_meta(L,reverse^=1);
_safe_meta(R,reverse^=1);}
if(to_push.a!=1||to_push.b!=0){_safe_meta(L,add_push(to_push));
_safe_meta(R,add_push(to_push));
to_push={};}}};}
#line 1 "cp-algo/structures/treap.hpp"
#line 1 "cp-algo/random/rng.hpp"
#include <chrono>
#include <random>
namespace cp_algo::random{std::mt19937_64 gen(
std::chrono::steady_clock::now().time_since_epoch().count()
);
uint64_t rng(){return gen();}}
#line 5 "cp-algo/structures/treap.hpp"
#include <array>
namespace cp_algo::structures::treap{template<typename meta>
struct node{using treap=node*;
meta _meta;
int prior=(int)random::rng();
size_t size=1;
treap children[2]={nullptr,nullptr};
enum subtree{L,R};
node(){}
node(meta _meta):_meta(_meta){}
node(meta _meta,int prior):_meta(_meta),prior(prior){}
static treap make_treap(auto...args){return new node(args...);}
treap pull(){_meta.pull(children[L],children[R]);
size=1+_safe(children[L],size)+_safe(children[R],size);
return this;}
treap push(){_meta.push(children[L],children[R]);
return this;}
treap set(subtree i,treap t){children[i]=t;
return pull();}
treap cut(subtree i){return children[i];}
static treap merge(treap A,treap B){if(!_safe(A,push())||!_safe(B,push())){return A?A:B;}else if(A->prior<B->prior){return A->set(R,merge(A->cut(R),B));}else{return B->set(L,merge(A,B->cut(L)));}}
static std::array<treap,2>split(treap A,size_t k){if(!_safe(A,push())){return{nullptr,nullptr};}else if(_safe(A->children[L],size)>=k){auto[split_L,split_R]=split(A->cut(L),k);
return{split_L,A->set(L,split_R)};}else{k-=_safe(A->children[L],size)+1;
auto[split_L,split_R]=split(A->cut(R),k);
return{A->set(R,split_L),split_R};}}
static void exec_on_segment(treap&A,size_t l,size_t r,auto func){auto[LM,R]=split(A,r);
auto[L,M]=split(LM,l);
func(M);
A=merge(L,merge(M,R));}
static void insert(treap&A,size_t pos,treap t){auto[L,R]=split(A,pos);
A=merge(L,merge(t,R));}
static void erase(treap&A,size_t pos){auto[L,MR]=split(A,pos);
auto[M,R]=split(MR,1);
delete M;
A=merge(L,R);}
static void exec_on_each(treap&A,auto func){if(A){exec_on_each(A->children[L],func);
func(A);
exec_on_each(A->children[R],func);}}
treap pull_all(){_safe(children[L],pull_all());
_safe(children[R],pull_all());
return pull();}
treap push_all(){push();
_safe(children[L],push_all());
_safe(children[R],push_all());
return this;}
static treap build(auto const&nodes){std::vector<treap>st;
for(auto cur:nodes){while(st.size()>=2&&st[st.size()-2]->prior>cur->prior){st.pop_back();}
if(!st.empty()&&st.back()->prior>cur->prior){cur->set(L,st.back());
st.pop_back();}
if(!st.empty()&&st.back()->prior<cur->prior){st.back()->set(R,cur);}
st.push_back(cur);}
return st.empty()?nullptr:st[0]->pull_all();}};
struct null_meta{void pull(auto const,auto const){}
void push(auto&,auto&){}};}
#line 6 "verify/structures/treap/dynamic_sequence_range_affine_range_sum.test.cpp"
#include <bits/stdc++.h>
using namespace std;
using namespace cp_algo::structures::treap;
using base=cp_algo::math::modint<998244353>;
using meta=metas::reverse_meta<base>;
using node_t=node<meta>;
using treap=node_t::treap;
void solve(){istream_iterator<int>input(cin);
int n=*input++;
int q=*input++;
vector<treap>nodes(n);
generate_n(begin(nodes),n,[&](){return node_t::make_treap(meta(*input++));});
auto me=node_t::build(nodes);
while(q--){int t=*input++;
if(t==0){int i=*input++;
base x=*input++;
node_t::insert(me,i,node_t::make_treap(meta(x)));}else if(t==1){node_t::erase(me,*input++);}else if(t==2){int l=*input++;
int r=*input++;
node_t::exec_on_segment(me,l,r,[](auto&t){_safe_meta(t,reverse=1);});}else if(t==3){int l=*input++;
int r=*input++;
base b=*input++;
base c=*input++;
node_t::exec_on_segment(me,l,r,[b,c](auto&t){_safe_meta(t,add_push(meta::lin(b,c)));});}else{int l=*input++;
int r=*input++;
node_t::exec_on_segment(me,l,r,[](auto t){cout<<_safe_meta(t,sum)<<"\n";});}}}
signed main(){ios::sync_with_stdio(0);
cin.tie(0);
int t=1;
while(t--){solve();}}