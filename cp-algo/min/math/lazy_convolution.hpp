#ifndef CP_ALGO_MATH_LAZY_MULTIPLY_HPP
#define CP_ALGO_MATH_LAZY_MULTIPLY_HPP
#include "convolution.hpp"
#include <algorithm>
namespace cp_algo::math{template<typename base>auto lazy_multiply(base a0,base b0,auto&&get_ab,size_t n){big_vector<base>A={a0},B={b0};big_vector<base>C(n);C[0]=a0*b0;auto cdq=[&](this auto&&cdq,size_t l,size_t r)->void{if(r-l==1){auto[al,bl]=get_ab(A,B,C,l);A.push_back(al);B.push_back(bl);C[l]+=A[l]*B[0]+A[0]*B[l];return;}auto m=(l+r)/2;cdq(l,m);auto A_pref=std::span(A).subspan(0,std::min(m,r-l));auto B_pref=std::span(B).subspan(0,std::min(l,r-l));big_vector<base>A_suf(std::from_range,std::span(A).subspan(l,m-l));big_vector<base>B_suf(std::from_range,std::span(B).subspan(l,m-l));convolution_prefix(A_suf,B_pref,r-l);convolution_prefix(B_suf,A_pref,r-l);A_suf.resize(r-l);B_suf.resize(r-l);for(size_t i=m;i<r;i++){C[i]+=A_suf[i-l]+B_suf[i-l];}cdq(m,r);};cdq(1,n);return C;}}
#endif
