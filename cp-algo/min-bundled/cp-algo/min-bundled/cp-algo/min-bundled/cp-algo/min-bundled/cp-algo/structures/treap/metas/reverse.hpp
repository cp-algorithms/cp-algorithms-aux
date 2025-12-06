#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/structures/treap/metas/reverse.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/structures/treap/metas/reverse.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/structures/treap/metas/reverse.hpp"
#line 1 "cp-algo/structures/treap/metas/reverse.hpp"
#line 1 "cp-algo/structures/treap/metas/base.hpp"
#line 1 "cp-algo/structures/treap/common.hpp"
#define _safe(t, op) (t ? t->op : typename std::remove_reference_t<decltype(t->op)>())
#line 4 "cp-algo/structures/treap/metas/base.hpp"
#include <functional>
#include <algorithm>
#include <cstdint>
#define _safe_meta(i, op) _safe(i, _meta.op)
namespace cp_algo::structures::treap::metas{struct base_meta{void pull(auto const,auto const){}
void push(auto&,auto&){}};}
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