#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/structures/segtree/metas/affine.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/structures/segtree/metas/affine.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/structures/segtree/metas/affine.hpp"
#line 1 "cp-algo/structures/segtree/metas/affine.hpp"
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