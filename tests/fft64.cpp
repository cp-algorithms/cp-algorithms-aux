#include "cp-algo/math/fft64.hpp"
#include <cassert>
#include <random>
using namespace cp_algo::math;
int main(){std::mt19937_64 rng(8147);for(size_t n:{1,2,7,8,9,31,32,33,64,65,129}){std::vector<uint64_t>a(n),b(n+1),want(2*n);for(auto&x:a)x=rng();for(auto&x:b)x=rng();for(size_t i=0;i<a.size();i++)for(size_t j=0;j<b.size();j++)want[i+j]+=a[i]*b[j];fft::conv64(a,b);assert(a==want);}std::vector<uint64_t>a{100000},b{70368744177665};fft::conv64(a,b);assert(a[0]==7036874417766500000ULL);}
