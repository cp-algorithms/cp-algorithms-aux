#ifndef CP_ALGO_STRUCTURES_SEGMENT_TREE_HPP
#define CP_ALGO_STRUCTURES_SEGMENT_TREE_HPP
#include <vector>
#include <numeric>
namespace cp_algo::structures{
template<typename meta>
struct segtree_t{
const size_t N;
std::vector<meta>_meta;
segtree_t(size_t n):N(n),_meta(4*N){}
segtree_t(std::vector<meta>leafs):N(size(leafs)),_meta(4*N){
build(leafs);
}
void pull(size_t v,size_t l,size_t r){
if(r-l>1){
_meta[v].pull(_meta[2*v],_meta[2*v+1],l,r);
}
}
void push(size_t v,size_t l,size_t r){
if(r-l>1){
_meta[v].push(&_meta[2*v],&_meta[2*v+1],l,r);
}else{
_meta[v].push(nullptr,nullptr,l,r);
}
}
void build(auto&a,size_t v,size_t l,size_t r){
if(r-l==1){
if(l<size(a)){
_meta[v]=a[l];
}
}else{
size_t m=std::midpoint(l,r);
build(a,2*v,l,m);
build(a,2*v+1,m,r);
pull(v,l,r);
}
}
void build(auto&a){
build(a,1,0,N);
}
void exec_on_segment(size_t a,size_t b,auto func,auto proceed,auto stop,size_t v,size_t l,size_t r){
push(v,l,r);
if(r<=a||b<=l||stop(_meta[v])){
return;
}else if(a<=l&&r<=b&&proceed(_meta[v])){
func(_meta[v]);
push(v,l,r);
}else{
size_t m=std::midpoint(l,r);
exec_on_segment(a,b,func,proceed,stop,2*v,l,m);
exec_on_segment(a,b,func,proceed,stop,2*v+1,m,r);
pull(v,l,r);
}
}
static constexpr auto default_true=[](auto const&){return true;};
static constexpr auto default_false=[](auto const&){return false;};
void exec_on_segment(size_t a,size_t b,auto func,auto proceed,auto stop){
exec_on_segment(a,b,func,proceed,stop,1,0,N);
}
void exec_on_segment(size_t a,size_t b,auto func){
exec_on_segment(a,b,func,default_true,default_false);
}
};
}
#endif
