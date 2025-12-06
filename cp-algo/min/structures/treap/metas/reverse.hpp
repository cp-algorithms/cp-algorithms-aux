#ifndef CP_ALGO_STRUCTURES_TREAP_METAS_REVERSE_HPP
#define CP_ALGO_STRUCTURES_TREAP_METAS_REVERSE_HPP
#include "base.hpp"
#include "../../../math/affine.hpp"
#include <algorithm>
namespace cp_algo::structures::treap::metas{
template<typename base>
struct reverse_meta:base_meta{
using lin=math::lin<base>;
base val;
size_t sz=1;
bool reverse=false;
base sum=val;
lin to_push={};
reverse_meta(base val):val(val){}
void pull(auto const L,auto const R){
sum=val+_safe_meta(L,sum)+_safe_meta(R,sum);
sz=1+_safe_meta(L,sz)+_safe_meta(R,sz);
}
void add_push(lin const&t){
val=t.eval(val);
sum=t.a*sum+t.b*sz;
to_push.prepend(t);
}
void push(auto&L,auto&R){
if(reverse){
reverse=false;
std::swap(L,R);
_safe_meta(L,reverse^=1);
_safe_meta(R,reverse^=1);
}
if(to_push.a!=1||to_push.b!=0){
_safe_meta(L,add_push(to_push));
_safe_meta(R,add_push(to_push));
to_push={};
}
}
};
}
#endif
