#ifndef CP_ALGO_MATH_MULTIVAR_INV_HPP
#define CP_ALGO_MATH_MULTIVAR_INV_HPP
#include "multivar.hpp"
#include <algorithm>
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::fft{template<modint_type base>multivar<base>multivar_inv(multivar<base>const&a,auto const&target_dim){assert(a.data[0]!=base(0));size_t K=a.dim.size();big_vector<size_t>init_dim(K,1);multivar<base>b(init_dim);b.data[0]=a.data[0].inv();size_t degree=1;size_t total_degree=0;for(size_t i=0;i<K;i++){total_degree+=target_dim[i];}big_vector<size_t>next_dim(K);while(degree<total_degree){degree*=2;for(size_t i=0;i<K;i++){next_dim[i]=std::min(target_dim[i],degree);}multivar<base>b_ext(next_dim);b_ext.assign_prefix_from(b);auto c=a.truncated(next_dim);c.mul(b_ext);for(auto&x:c.data){x=-x;}c.data[0]+=base(2);b_ext.mul(c);b=std::move(b_ext);}b.truncate_inplace(target_dim);return b;}}
#pragma GCC pop_options
#endif
