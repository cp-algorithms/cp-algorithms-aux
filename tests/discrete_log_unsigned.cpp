#include "cp-algo/number_theory/discrete_log.hpp"
#include <cassert>
#include <random>
using namespace cp_algo::math;
int main() {
    assert(discrete_log<uint64_t>(3,9,17)==2);
    assert(discrete_log<uint32_t>(2,10,15,5)==1);
    std::mt19937 rng(67281);
    for(uint64_t m=1;m<=100;m++) for(int trial=0;trial<100;trial++) {
        uint64_t a=rng()%m,b=rng()%m,c=rng()%m;
        std::optional<uint64_t> want;
        uint64_t value=a;
        for(uint64_t x=0;x<=m;x++) {if(value==c){want=x;break;}value=value*b%m;}
        assert(discrete_log(b,c,m,a)==want);
        assert(discrete_log(uint32_t(b),uint32_t(c),uint32_t(m),uint32_t(a))==want);
    }
}
