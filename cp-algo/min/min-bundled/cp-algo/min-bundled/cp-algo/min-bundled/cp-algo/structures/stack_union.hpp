#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/structures/stack_union.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/structures/stack_union.hpp"
#line 1 "cp-algo/structures/stack_union.hpp"
#line 1 "cp-algo/util/big_alloc.hpp"
#include <vector>
#include <cstddef>
#include <iostream>
#if defined(__linux__) || defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#  define CP_ALGO_USE_MMAP 1
#  include <sys/mman.h>
#else
#  define CP_ALGO_USE_MMAP 0
#endif
namespace cp_algo{template<typename T,std::size_t Align=32>
class big_alloc{static_assert(Align>=alignof(void*),"Align must be at least pointer-size");
static_assert(std::popcount(Align)==1,"Align must be a power of two");
public:using value_type=T;
template<class U>struct rebind{using other=big_alloc<U,Align>;};
constexpr bool operator==(const big_alloc&)const=default;
constexpr bool operator!=(const big_alloc&)const=default;
big_alloc()noexcept=default;
template<typename U,std::size_t A>
big_alloc(const big_alloc<U,A>&)noexcept{}
[[nodiscard]]T*allocate(std::size_t n){std::size_t padded=round_up(n*sizeof(T));
std::size_t align=std::max<std::size_t>(alignof(T),Align);
#if CP_ALGO_USE_MMAP
if(padded>=MEGABYTE){void*raw=mmap(nullptr,padded,
PROT_READ|PROT_WRITE,
MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
madvise(raw,padded,MADV_HUGEPAGE);
madvise(raw,padded,MADV_POPULATE_WRITE);
return static_cast<T*>(raw);}
#endif
return static_cast<T*>(::operator new(padded,std::align_val_t(align)));}
void deallocate(T*p,std::size_t n)noexcept{if(!p)return;
std::size_t padded=round_up(n*sizeof(T));
std::size_t align=std::max<std::size_t>(alignof(T),Align);
#if CP_ALGO_USE_MMAP
if(padded>=MEGABYTE){munmap(p,padded);return;}
#endif::operator delete(p,padded,std::align_val_t(align));}
private:static constexpr std::size_t MEGABYTE=1<<20;
static constexpr std::size_t round_up(std::size_t x)noexcept{return(x+Align-1)/Align*Align;}};
template<typename T>
using big_vector=std::vector<T,big_alloc<T>>;}
#line 5 "cp-algo/structures/stack_union.hpp"
#include <iterator>
#include <ranges>
namespace cp_algo::structures{template<class datatype>
struct stack_union{stack_union(int n=0):head(n),next(1),data(1){}
void push(int v,datatype const&vdata){next.push_back(head[v]);
head[v]=(int)std::size(next)-1;
data.push_back(vdata);}
template<typename... Args>
void emplace(int v,Args&&... vdata){next.push_back(head[v]);
head[v]=(int)std::size(next)-1;
data.emplace_back(std::forward<Args>(vdata)...);}
void reserve(int m){data.reserve(m);
next.reserve(m);}
size_t size()const{return std::size(head);}
size_t nodes()const{return std::size(data);}
template<typename Su>
struct _iterator{using value_type=std::conditional_t<std::is_const_v<Su>,const datatype,datatype>;
using difference_type=std::ptrdiff_t;
Su*su=nullptr;
int sv=0;
value_type&operator*()const{return su->data[sv];}
_iterator&operator++(){sv=su->next[sv];
return*this;}
_iterator operator++(int){auto tmp=*this;++*this;return tmp;}
friend bool operator==(_iterator const&it,std::default_sentinel_t){return it.sv==0;}};
using iterator=_iterator<stack_union<datatype>>;
using const_iterator=_iterator<const stack_union<datatype>>;
auto operator[](this auto&&self,int v){using Iter=_iterator<std::remove_reference_t<decltype(self)>>;
return std::ranges::subrange(Iter{&self,self.head[v]},std::default_sentinel);}
big_vector<int>head,next;
big_vector<datatype>data;};}