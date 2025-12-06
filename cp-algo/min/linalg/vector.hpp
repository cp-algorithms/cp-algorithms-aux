#ifndef CP_ALGO_LINALG_VECTOR_HPP
#define CP_ALGO_LINALG_VECTOR_HPP
#include "../random/rng.hpp"
#include "../number_theory/modint.hpp"
#include "../util/big_alloc.hpp"
#include "../util/simd.hpp"
#include "../util/checkpoint.hpp"
#include <functional>
#include <algorithm>
#include <valarray>
#include <iostream>
#include <iterator>
#include <cassert>
#include <ranges>
namespace cp_algo::linalg{template<typename base,class Alloc=big_alloc<base>>
struct vec:std::basic_string<base,std::char_traits<base>,Alloc>{using Base=std::basic_string<base,std::char_traits<base>,Alloc>;
using Base::Base;
vec(Base const&t):Base(t){}
vec(Base&&t):Base(std::move(t)){}
vec(size_t n):Base(n,base()){}
vec(auto&&r):Base(std::ranges::to<Base>(r)){}
static vec ei(size_t n,size_t i){vec res(n);
res[i]=1;
return res;}
auto operator-()const{return*this|std::views::transform([](auto x){return-x;});}
auto operator*(base t)const{return*this|std::views::transform([t](auto x){return x*t;});}
auto operator*=(base t){for(auto&it:*this){it*=t;}
return*this;}
virtual void add_scaled(vec const&b,base scale,size_t i=0){if(scale!=base(0)){for(;i<size(*this);i++){(*this)[i]+=scale*b[i];}}}
virtual vec const&normalize(){return static_cast<vec&>(*this);}
virtual base normalize(size_t i){return(*this)[i];}
void read(){for(auto&it:*this){std::cin>>it;}}
void print()const{for(auto&it:*this){std::cout<<it<<" ";}
std::cout<<"\n";}
static vec random(size_t n){vec res(n);
std::ranges::generate(res,random::rng);
return res;}
vec operator|(vec const&t)const{return std::views::join(std::array{std::views::all(*this),
std::views::all(t)});}
std::pair<size_t,base>find_pivot(){if(pivot==size_t(-1)){pivot=0;
while(pivot<size(*this)&&normalize(pivot)==base(0)){pivot++;}
if(pivot<size(*this)){pivot_inv=base(1)/(*this)[pivot];}}
return{pivot,pivot_inv};}
void reduce_by(vec&t){auto[pivot,pinv]=t.find_pivot();
if(pivot<size(*this)){add_scaled(t,-normalize(pivot)*pinv,pivot);}}
private:size_t pivot=-1;
base pivot_inv;};
template<math::modint_type base,class Alloc=big_alloc<base>>
struct modint_vec:vec<base,Alloc>{using Base=vec<base,Alloc>;
using Base::Base;
modint_vec(Base const&t):Base(t){}
modint_vec(Base&&t):Base(std::move(t)){}
void add_scaled(Base const&b,base scale,size_t i=0)override{static_assert(base::bits>=64,"Only wide modint types for linalg");
if(scale!=base(0)){assert(Base::size()==b.size());
size_t n=size(*this);
u64x4 scaler=u64x4()+scale.getr();
if(is_aligned(&(*this)[0])&&is_aligned(&b[0]))
for(i-=i%4;i+3<n;i+=4){auto&ai=vector_cast<u64x4>((*this)[i]);
auto bi=vector_cast<u64x4 const>(b[i]);
#ifdef __AVX2__
ai+=u64x4(_mm256_mul_epu32(__m256i(scaler),__m256i(bi)));
#else
ai+=scaler*bi;}
for(;i<n;i++){(*this)[i].add_unsafe(b[i].getr_direct()*scale.getr());}
if(++counter==4){for(auto&it:*this){it.pseudonormalize();}
counter=0;}}}
Base const&normalize()override{for(auto&it:*this){it.normalize();}
return*this;}
base normalize(size_t i)override{return(*this)[i].normalize();}
private:size_t counter=0;};}
#endif
#endif