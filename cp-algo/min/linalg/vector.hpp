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
#include <array>
#include <cstring>
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::linalg{template<typename base,class Alloc=big_alloc<base>>struct vec:std::basic_string<base,std::char_traits<base>,Alloc>{using Base=std::basic_string<base,std::char_traits<base>,Alloc>;using Base::Base;vec(Base const&t):Base(t){}vec(Base&&t):Base(std::move(t)){}vec(size_t n):Base(n,base()){}vec(auto&&r):Base(r.begin(),r.end()){}static vec ei(size_t n,size_t i){vec res(n);res[i]=1;return res;}auto operator-()const{return*this|std::views::transform([](auto x){return-x;});}auto operator*(base t)const{return*this|std::views::transform([t](auto x){return x*t;});}vec&operator*=(base t){for(auto&it:*this){it*=t;}return*this;}virtual void add_scaled(vec const&b,base scale,size_t i=0){if(scale!=base(0)){for(;i<size(*this);i++){(*this)[i]+=scale*b[i];}}}virtual vec const&normalize(){return static_cast<vec&>(*this);}virtual base normalize(size_t i){return(*this)[i];}void read(){for(auto&it:*this){std::cin>>it;}}void print()const{for(auto&it:*this){std::cout<<it<<' ';}std::cout<<'\n';}static vec random(size_t n){vec res(n);std::ranges::generate(res,random::rng);return res;}vec operator|(vec const&t)const{return std::views::join(std::array{std::views::all(*this),std::views::all(t)});}std::pair<size_t,base>find_pivot(){if(pivot==size_t(-1)){pivot=0;while(pivot<size(*this)&&normalize(pivot)==base(0)){pivot++;}if(pivot<size(*this)){pivot_inv=base(1)/(*this)[pivot];}}return{pivot,pivot_inv};}void reduce_by(vec&t){auto[pivot,pinv]=t.find_pivot();if(pivot<size(*this)){add_scaled(t,-normalize(pivot)*pinv,pivot);}}private:size_t pivot=-1;base pivot_inv;};template<math::modint_type base,class Alloc=big_alloc<base>>struct modint_vec:vec<base,Alloc>{using Base=vec<base,Alloc>;using Base::Base;modint_vec(Base const&t):Base(t){}modint_vec(Base&&t):Base(std::move(t)){}void add_scaled(Base const&b,base scale,size_t i=0)override{static_assert(base::bits>=64,"Only wide modint types for linalg");if(scale!=base(0)){assert(Base::size()==b.size());size_t n=size(*this);u64x4 scaler=u64x4()+scale.getr();bool aligned=is_aligned(&(*this)[0])&&is_aligned(&b[0]);if(aligned)i-=i%4;bool reduce=++counter==accumulation_period();if(reduce){counter=0;for(size_t j=0;j<i;j++)(*this)[j].pseudonormalize();}if(aligned)for(;i+3<n;i+=4){auto&ai=vector_cast<u64x4>((*this)[i]);auto bi=vector_cast<u64x4 const>(b[i]);
#ifdef __AVX2__
ai+=u64x4(_mm256_mul_epu32(__m256i(scaler),__m256i(bi)));
#else
ai+=scaler*bi;
#endif
if(reduce)ai=shrink(ai);}for(;i<n;i++){(*this)[i].add_unsafe(b[i].getr_direct()*scale.getr());if(reduce)(*this)[i].pseudonormalize();}}}Base const&normalize()override{for(auto&it:*this){it.normalize();}return*this;}base normalize(size_t i)override{return(*this)[i].normalize();}private:template<typename,typename>friend struct matrix;static size_t accumulation_period(){return base::remod()==base::mod()&&base::mod()<(1LL<<30)?8:4;}static u64x4 mul(u64x4 a,u64x4 b){
#ifdef __AVX2__
return u64x4(_mm256_mul_epu32(__m256i(a),__m256i(b)));
#else
return a*b;
#endif
}static u64x4 shrink(u64x4 a){auto b=a-(u64x4()+base::modmod8());return a<b?a:b;}static void add_scaled_pair(modint_vec&x,modint_vec&y,Base const&p,Base const&q,std::array<base,4>c,size_t first=0){if(std::ranges::find(c,base(0))!=c.end()){x.add_scaled(p,c[0],first);x.add_scaled(q,c[1],first);y.add_scaled(p,c[2],first);y.add_scaled(q,c[3],first);return;}size_t n=x.size();assert(y.size()==n&&p.size()==n&&q.size()==n&&first<=n);size_t period=accumulation_period();auto prepare=[&](modint_vec&a){if(a.counter+2>period){for(auto&v:a)v.pseudonormalize();a.counter=0;}a.counter+=2;if(a.counter!=period)return false;a.counter=0;for(size_t i=0;i<first;i++)a[i].pseudonormalize();return true;};bool nx=prepare(x),ny=prepare(y);auto*__restrict__ dx=x.data();auto*__restrict__ dy=y.data();auto const*__restrict__ sp=p.data();auto const*__restrict__ sq=q.data();uint64_t xp=c[0].getr(),xq=c[1].getr(),yp=c[2].getr(),yq=c[3].getr();u64x4 xp4=u64x4()+xp,xq4=u64x4()+xq;u64x4 yp4=u64x4()+yp,yq4=u64x4()+yq;size_t i=first;for(;i+4<=n;i+=4){u64x4 vx,vy,vp,vq;std::memcpy(&vx,dx+i,sizeof vx);std::memcpy(&vy,dy+i,sizeof vy);std::memcpy(&vp,sp+i,sizeof vp);std::memcpy(&vq,sq+i,sizeof vq);vx+=mul(xp4,vp)+mul(xq4,vq);vy+=mul(yp4,vp)+mul(yq4,vq);if(nx)vx=shrink(vx);if(ny)vy=shrink(vy);std::memcpy(dx+i,&vx,sizeof vx);std::memcpy(dy+i,&vy,sizeof vy);}for(;i<n;i++){dx[i].add_unsafe(xp*sp[i].getr_direct()+xq*sq[i].getr_direct());dy[i].add_unsafe(yp*sp[i].getr_direct()+yq*sq[i].getr_direct());if(nx)dx[i].pseudonormalize();if(ny)dy[i].pseudonormalize();}}size_t counter=0;};template<typename base>requires(base::bits<=32)struct modint_vec<base>:vec<base>{using Base=vec<base>;using Base::Base;modint_vec(Base const&t):Base(t){}modint_vec(Base&&t):Base(std::move(t)){}void add_scaled(Base const&b,base scale,size_t first=0)override{if(scale==base(0))return;if(&b==this)Base::add_scaled(b,scale,first);else add_scaled_batch<1,1>({this},{&b},{scale},first);}private:template<typename,typename>friend struct matrix;static constexpr size_t batch_size=8;static constexpr bool use_simd=[]{if constexpr(requires{std::integral_constant<uint32_t,base::mod()>{};std::integral_constant<uint32_t,base::remod()>{};}){return sizeof(base)==sizeof(uint32_t)&&base::mod()>1&&base::mod()%2&&base::mod()<(1U<<30)&&base::remod()==base::mod();}else return false;}();static u64x4 mul(u64x4 a,u64x4 b){
#ifdef __AVX2__
return u64x4(_mm256_mul_epu32(__m256i(a),__m256i(b)));
#else
return low32(a)*low32(b);
#endif
}template<size_t count,size_t rows>static void add_scaled_batch(std::array<modint_vec*,rows>const&dst,std::array<Base const*,count>const&src,std::array<base,rows*count>const&c,size_t first=0){static_assert(count<=batch_size&&(rows==1||rows==2));if constexpr(use_simd){constexpr uint32_t mod=base::mod(),inv=math::inv2(uint32_t(-mod));std::array<uint32_t,rows*count>scale;for(size_t t=0;t<rows*count;t++){scale[t]=uint32_t((uint64_t(c[t].getr())<<32)%mod);}auto*__restrict__ dx=dst[0]->data();auto*__restrict__ dy=rows==2?dst[1]->data():nullptr;size_t n=dst[0]->size();for(;first+8<=n;first+=8){u64x4 acc[rows][2]{};
#pragma GCC unroll 1
for(size_t t=0;t<count;t++){u64x4 p;std::memcpy(&p,src[t]->data()+first,sizeof p);auto q=p>>32;for(size_t row=0;row<rows;row++){auto v=u64x4(u32x8()+scale[row*count+t]);acc[row][0]+=mul(p,v);acc[row][1]+=mul(q,v);}}for(size_t row=0;row<rows;row++){for(auto&v:acc[row])v=montgomery_reduce(v,mod,inv);auto*out=(row?dy:dx)+first;u32x8 old;std::memcpy(&old,out,sizeof old);auto z=old+u32x8(acc[row][0]|(acc[row][1]<<32));z=z<z-2*mod?z:z-2*mod;z=z<z-mod?z:z-mod;std::memcpy(out,&z,sizeof z);}}}for(size_t t=0;t<count;t++)for(size_t row=0;row<rows;row++){dst[row]->Base::add_scaled(*src[t],c[row*count+t],first);}}static void add_scaled_pair(modint_vec&x,modint_vec&y,Base const&p,Base const&q,std::array<base,4>c,size_t first=0){add_scaled_batch<2,2>({&x,&y},{&p,&q},c,first);}};}
#pragma GCC pop_options
#endif
