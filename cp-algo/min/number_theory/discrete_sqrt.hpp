#ifndef CP_ALGO_NUMBER_THEORY_DISCRETE_SQRT_HPP
#define CP_ALGO_NUMBER_THEORY_DISCRETE_SQRT_HPP
#include "modint.hpp"
#include "../random/rng.hpp"
#include "../math/affine.hpp"
namespace cp_algo::math{
template<modint_type base>
std::optional<base>sqrt(base b){
if(b==base(0)){
return base(0);
}else if(bpow(b,(b.mod()-1)/2)!=base(1)){
return std::nullopt;
}else{
while(true){
base z=random::rng();
if(z*z==b){
return z;
}
lin<base>x(1,z,b);
x=bpow(x,(b.mod()-1)/2,lin<base>(0,1,b));
if(x.a!=base(0)){
return x.a.inv();
}
}
}
}
}
#endif
