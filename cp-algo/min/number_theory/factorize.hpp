#ifndef CP_ALGO_MATH_FACTORIZE_HPP
#define CP_ALGO_MATH_FACTORIZE_HPP
#include "primality.hpp"
#include "../util/generator.hpp"
#include "../random/rng.hpp"
namespace cp_algo::math{template<typename _Int>auto proper_divisor(_Int m){using Int=std::make_signed_t<_Int>;using base=dynamic_modint<Int>;return m%2==0?2:base::with_mod(m,[&](){base t=random::rng();auto f=[&](auto x){return x*x+t;};const base zero=0,one=1;base x=zero,y=zero;base g=one;while(g==one){for(int i=0;i<64;i++){x=f(x);y=f(f(y));if(x==y)[[unlikely]]{t=random::rng();x=y=zero;}else{base t=g*(x-y);g=t==zero?g:t;}}g=std::gcd(g.getr(),m);}return g.getr();});}template<typename Int>big_generator<Int>factorize(Int m){if(is_prime(m)){co_yield m;}else if(m>1){auto g=proper_divisor(m);co_yield std::ranges::elements_of(factorize(g));co_yield std::ranges::elements_of(factorize(m/g));}}template<typename Int>big_generator<Int>divisors_sqrt(Int m){for(Int i=1;i*i<=m;i++){if(m%i==0){co_yield i;if(i*i!=m){co_yield m/i;}}}}}
#endif
