#line 1 "verify/structures/segtree/range_chmin_chmax_add_range_sum.test.cpp"
#define PROBLEM "https:#line 1 "cp-algo/structures/segtree/metas/chmin_chmax_add.hpp"
#line 1 "cp-algo/structures/segtree/metas/base.hpp"
#include <cstddef>
namespace cp_algo::structures::segtree::metas{template<typename derived_meta>
struct base_meta{using meta=derived_meta;
virtual void pull(meta const&,meta const&,size_t,size_t){};
virtual void push(meta*,meta*,size_t,size_t){};};}
#line 4 "cp-algo/structures/segtree/metas/chmin_chmax_add.hpp"
#include <functional>
#include <algorithm>
#include <cstdint>
namespace cp_algo::structures::segtree::metas{struct chmin_chmax_sum_meta:base_meta<chmin_chmax_sum_meta>{static constexpr int64_t inf=1e12;
using meta=chmin_chmax_sum_meta;
int64_t sum=0,add=0;
template<typename Comp>
struct data{int64_t val;
size_t count=1;
int64_t second=std::max(inf,-inf,comp);
static const Comp comp;
data combine(data const&t)const{return comp(val,t.val)?data{val,count,std::min(second,t.val,comp)}:comp(t.val,val)?data{t.val,t.count,std::min(t.second,val,comp)}:data{val,count+t.count,std::min(second,t.second,comp)};}
void add(int64_t b){val+=b;
second+=b;}
int64_t normalize(int64_t L,int64_t R){int64_t old_val=val;
val=std::clamp(val,L,R);
second=std::clamp(second,L,R);
return count*(val-old_val);}
bool stop(int64_t b)const{return!comp(val,b);}
bool proceed(int64_t b)const{return comp(b,second);}};
data<std::less<>>mn={sum};
data<std::greater<>>mx={sum};
int64_t chmin=inf,chmax=-inf;
chmin_chmax_sum_meta(){}
chmin_chmax_sum_meta(int64_t val):sum(val){}
void pull(meta const&L,meta const&R,size_t,size_t)override{sum=L.sum+R.sum;
mn=L.mn.combine(R.mn);
mx=L.mx.combine(R.mx);}
void push(meta&t){t.add+=add;t.chmin+=add;t.chmax+=add;
t.chmin=std::clamp(t.chmin,chmax,chmin);
t.chmax=std::clamp(t.chmax,chmax,chmin);}
void push(meta*L,meta*R,size_t l,size_t r)override{if(r-l>1){push(*L);
push(*R);}
if(add){sum+=(r-l)*add;
mn.add(add);
mx.add(add);}
bool same=mn.val==mx.val;
auto to_add=mn.normalize(chmax,chmin)+mx.normalize(chmax,chmin);
sum+=same?to_add/2:to_add;
if(mn.val==mx.val){mx={mx.val,r-l};
mn={mn.val,r-l};}
add=0;
chmin=inf;
chmax=-inf;}
static auto proceed_chmin(int64_t b){return[b](meta const&t){return t.mx.proceed(b);};}
static auto stop_chmin(int64_t b){return[b](meta const&t){return t.mx.stop(b);};}
static auto proceed_chmax(int64_t b){return[b](meta const&t){return t.mn.proceed(b);};}
static auto stop_chmax(int64_t b){return[b](meta const&t){return t.mn.stop(b);};}};}
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
#line 5 "verify/structures/segtree/range_chmin_chmax_add_range_sum.test.cpp"
#include <bits/stdc++.h>
using namespace std;
using namespace cp_algo::structures;
using meta=segtree::metas::chmin_chmax_sum_meta;
void solve(){int n,q;
cin>>n>>q;
vector<meta>a(n);
for(int i=0;i<n;i++){int64_t ai;
cin>>ai;
a[i]={ai};}
segtree_t<meta>me(a);
while(q--){int t,l,r;
int64_t b;
cin>>t>>l>>r;
if(t==0){cin>>b;
me.exec_on_segment(l,r,
[b](auto&meta){meta.chmin=b;},
meta::proceed_chmin(b),meta::stop_chmin(b));}else if(t==1){cin>>b;
me.exec_on_segment(l,r,
[b](auto&meta){meta.chmax=b;},
meta::proceed_chmax(b),meta::stop_chmax(b));}else if(t==2){cin>>b;
me.exec_on_segment(l,r,
[b](auto&meta){meta.add=b;});}else{int64_t ans=0;
me.exec_on_segment(l,r,[&](auto&meta){ans+=meta.sum;});
cout<<ans<<"\n";}}}
signed main(){ios::sync_with_stdio(0);
cin.tie(0);
int t=1;
while(t--){solve();}}