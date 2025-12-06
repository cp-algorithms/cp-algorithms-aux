#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/min/structures/segtree/metas/chmin_chmax_add.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/min/structures/segtree/metas/chmin_chmax_add.hpp"
#line 1 "cp-algo/min/structures/segtree/metas/chmin_chmax_add.hpp"
#line 1 "cp-algo/min/structures/segtree/metas/base.hpp"
#include <cstddef>
namespace cp_algo::structures::segtree::metas{template<typename derived_meta>
struct base_meta{using meta=derived_meta;
virtual void pull(meta const&,meta const&,size_t,size_t){};
virtual void push(meta*,meta*,size_t,size_t){};};}
#line 4 "cp-algo/min/structures/segtree/metas/chmin_chmax_add.hpp"
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