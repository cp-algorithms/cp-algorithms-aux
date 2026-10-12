#include "cp-algo/number_theory/factorize.hpp"
#include <cassert>
#include <set>
using namespace cp_algo::math;
int main(){for(int n:{1,2,4,36,46340*46340,2147483647}){std::set<int> got,want;for(auto d:divisors_sqrt(n))got.insert(d);for(int64_t i=1;i*i<=n;i++)if(n%i==0){want.insert(i);want.insert(n/i);}assert(got==want);}}
