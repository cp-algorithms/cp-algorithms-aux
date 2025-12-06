#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/min/structures/csr.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/min/structures/csr.hpp"
#line 1 "cp-algo/min/structures/csr.hpp"
#line 1 "cp-algo/min/util/big_alloc.hpp"
#include <vector>
#include <cstddef>
#include <iostream>
#if defined(__linux__) || defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#  define CP_ALGO_USE_MMAP 1
#  include <sys/mman.h>
#else
#  define CP_ALGO_USE_MMAP 0
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
return static_cast<T*>(::operator new(padded,std::align_val_t(align)));}
void deallocate(T*p,std::size_t n)noexcept{if(!p)return;
std::size_t padded=round_up(n*sizeof(T));
std::size_t align=std::max<std::size_t>(alignof(T),Align);
#if CP_ALGO_USE_MMAP
if(padded>=MEGABYTE){munmap(p,padded);return;}::operator delete(p,padded,std::align_val_t(align));}
private:static constexpr std::size_t MEGABYTE=1<<20;
static constexpr std::size_t round_up(std::size_t x)noexcept{return(x+Align-1)/Align*Align;}};
template<typename T>
using big_vector=std::vector<T,big_alloc<T>>;}
#endif
#endif
#endif
#line 4 "cp-algo/min/structures/csr.hpp"
#include <algorithm>
#include <numeric>
#include <ranges>
#include <span>
namespace cp_algo::structures{template<class datatype>
struct csr{csr():head{0}{}
void new_row(){head.push_back(head.back());}
void push(datatype const&value){data.push_back(value);
head.back()++;}
template<typename... Args>
void emplace(Args&&... args){data.emplace_back(std::forward<Args>(args)...);
head.back()++;}
void reserve_rows(int n){head.reserve(n+1);}
void reserve_data(size_t n){data.reserve(n);}
void reverse_rows(){std::ranges::reverse(data);
std::adjacent_difference(head.begin(),head.end(),head.begin());
std::ranges::reverse(head|std::views::drop(1));
std::partial_sum(head.begin(),head.end(),head.begin());}
size_t size()const{return head.size()-1;}
size_t total_size()const{return data.size();}
auto operator[](this auto&&self,auto row){return std::span(self.data).subspan(self.head[row],self.head[row+1]-self.head[row]);}
auto rows(this auto&&self){return std::views::iota(size_t(0),self.size())
|std::views::transform([&self](auto i){return self[i];});}
big_vector<int>head;
big_vector<datatype>data;};}