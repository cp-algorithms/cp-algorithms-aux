#ifndef CP_ALGO_UTIL_BUMP_ALLOC_HPP
#define CP_ALGO_UTIL_BUMP_ALLOC_HPP
#include <cstddef>
#include "big_alloc.hpp"
namespace cp_algo{template<class T,size_t max_len>
struct bump_alloc{static char*buf;
static size_t buf_ind;
using value_type=T;
template<class U>struct rebind{using other=bump_alloc<U,max_len>;};
constexpr bool operator==(const bump_alloc&)const=default;
constexpr bool operator!=(const bump_alloc&)const=default;
bump_alloc()=default;
template<class U>bump_alloc(const U&){}
T*allocate(size_t n){buf_ind-=n*sizeof(T);
buf_ind&=0-alignof(T);
return(T*)(buf+buf_ind);}
void deallocate(T*,size_t){}};
template<class T,size_t max_len>
char*bump_alloc<T,max_len>::buf=big_alloc<char>().allocate(max_len*sizeof(T));
template<class T,size_t max_len>
size_t bump_alloc<T,max_len>::buf_ind=max_len*sizeof(T);}
#endif