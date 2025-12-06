#ifndef CP_ALGO_STRUCTURES_STACK_UNION_HPP
#define CP_ALGO_STRUCTURES_STACK_UNION_HPP
#include "../util/big_alloc.hpp"
#include <cstddef>
#include <iterator>
#include <ranges>
namespace cp_algo::structures{
template<class datatype>
struct stack_union{
stack_union(int n=0):head(n),next(1),data(1){}
void push(int v,datatype const&vdata){
next.push_back(head[v]);
head[v]=(int)std::size(next)-1;
data.push_back(vdata);
}
template<typename... Args>
void emplace(int v,Args&&... vdata){
next.push_back(head[v]);
head[v]=(int)std::size(next)-1;
data.emplace_back(std::forward<Args>(vdata)...);
}
void reserve(int m){
data.reserve(m);
next.reserve(m);
}
size_t size()const{return std::size(head);}
size_t nodes()const{return std::size(data);}
template<typename Su>
struct _iterator{
using value_type=std::conditional_t<std::is_const_v<Su>,const datatype,datatype>;
using difference_type=std::ptrdiff_t;
Su*su=nullptr;
int sv=0;
value_type&operator*()const{return su->data[sv];}
_iterator&operator++(){
sv=su->next[sv];
return*this;
}
_iterator operator++(int){auto tmp=*this;++*this;return tmp;}
friend bool operator==(_iterator const&it,std::default_sentinel_t){
return it.sv==0;
}
};
using iterator=_iterator<stack_union<datatype>>;
using const_iterator=_iterator<const stack_union<datatype>>;
auto operator[](this auto&&self,int v){
using Iter=_iterator<std::remove_reference_t<decltype(self)>>;
return std::ranges::subrange(Iter{&self,self.head[v]},std::default_sentinel);
}
big_vector<int>head,next;
big_vector<datatype>data;
};
}
#endif
