#ifndef CP_ALGO_MATH_FFT_SIMPLE_HPP
#define CP_ALGO_MATH_FFT_SIMPLE_HPP
#include "../random/rng.hpp"
#include "../math/common.hpp"
#include "../math/cvector.hpp"
CP_ALGO_SIMD_PRAGMA_PUSHnamespace cp_algo::math::fft{struct dft_simple{cp_algo::math::fft::cvector cv;dft_simple(auto const&a,size_t n):cv(n){for(size_t i=0;i<std::min(std::size(a),n);i++){real(cv.at(i))[i%4]=ftype(a[i]);imag(cv.at(i))[i%4]=ftype(i+n<std::size(a)?a[i+n]:0);}checkpoint("dft64 init");cv.fft();}void dot(dft_simple const&t){cv.dot(t.cv);}void recover_mod(auto&res,size_t k){cv.ifft();size_t n=cv.size();for(size_t i=0;i<std::min(k,n);i++){res[i]=llround(real(cv.get(i)));}for(size_t i=n;i<k;i++){res[i]=llround(imag(cv.get(i-n)));}cp_algo::checkpoint("recover mod");}};void conv_simple(auto&a,auto const&b){if(empty(a)||empty(b)){a.clear();return;}size_t n=a.size(),m=b.size();size_t N=std::max(flen,std::bit_ceil(n+m-1)/2);dft_simple A(a,N),B(b,N);A.dot(B);a.resize(n+m-1);A.recover_mod(a,n+m-1);}}
#pragma GCC pop_options
#endif
