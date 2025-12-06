#ifndef CP_ALGO_NUMBER_THEORY_EULER_HPP
#define CP_ALGO_NUMBER_THEORY_EULER_HPP
#include "factorize.hpp"
namespace cp_algo::math{auto euler_phi(auto m){auto primes=to<std::vector>(factorize(m));
std::ranges::sort(primes);
auto[from,to]=std::ranges::unique(primes);
primes.erase(from,to);
auto ans=m;
for(auto it:primes){ans-=ans/it;}
return ans;}
template<modint_type base>
auto period(base x){auto ans=euler_phi(base::mod());
base x0=bpow(x,ans);
for(auto t:factorize(ans)){while(ans%t==0&&x0*bpow(x,ans/t)==x0){ans/=t;}}
return ans;}
template<typename _Int>
_Int primitive_root(_Int p){using Int=std::make_signed_t<_Int>;
using base=dynamic_modint<Int>;
return base::with_mod(p,[p](){base t=1;
while(period(t)!=p-1){t=random::rng();}
return t.getr();});}}
#endif