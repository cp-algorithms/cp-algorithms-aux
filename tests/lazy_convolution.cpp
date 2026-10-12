#include <bits/stdc++.h>
#include "cp-algo/math/combinatorics.hpp"
#include "cp-algo/math/lazy_convolution.hpp"
#include "cp-algo/math/poly/series/exp.hpp"
using namespace cp_algo;using namespace cp_algo::math;using T=modint<998244353>;
int main(){std::cout<<std::unitbuf;std::mt19937 rng(6517);int checks=0;for(int n:{2,3,4,5,7,8,9,15,16,17,31,32,33,63,64,65,127,128,129,257})for(int rep=0;rep<4;rep++){
 big_vector<T>a(n),b(n),want(n);for(auto&x:a)x=rng()%T::mod();for(auto&x:b)x=rng()%T::mod();for(int i=0;i<n;i++)for(int j=0;j+i<n;j++)want[i+j]+=a[i]*b[j];auto c=lazy_multiply<T>(a[0],b[0],[&](auto&,auto&,auto&,size_t k){return std::pair{a[k],b[k]};},n);if(c!=want)throw std::runtime_error("full static");checks++;
 c=lazy_multiply<T>(a,b[0],[&](auto&,auto&,size_t k){return b[k];},n);if(c!=want)throw std::runtime_error("semi static");checks++;
 a[0]=2;b[0]=3;want[0]=6;for(int k=1;k<n;k++){a[k]=want[k-1]+T(k);b[k]=want[k-1]-T(k);want[k]=0;for(int i=0;i<=k;i++)want[k]+=a[i]*b[k-i];}
 c=lazy_multiply<T>(T(2),T(3),[](auto&,auto&,auto&prev,size_t k){return std::pair{prev[k-1]+T(k),prev[k-1]-T(k)};},n);if(c!=want)throw std::runtime_error("full feedback");checks++;
 b[0]=3;want[0]=a[0]*b[0];for(int k=1;k<n;k++){b[k]=want[k-1];want[k]=0;for(int i=0;i<=k;i++)want[k]+=a[i]*b[k-i];}c=lazy_multiply<T>(a,T(3),[](auto&,auto&prev,size_t k){return prev[k-1];},n);if(c!=want)throw std::runtime_error("semi feedback");checks++;
 a[0]=0;auto e=lazy_exp<T>(a[1],[&](auto&,auto&,size_t k){return k<a.size()?a[k]:T(0);},n);auto expected=exp(poly_t<T>(a),n);for(int k=0;k<n;k++)if(e[k]!=expected[k])throw std::runtime_error("exp");checks++;
}std::cout<<"PASS "<<checks<<" lazy callback/feedback checks\n";    for(size_t n:{size_t(0),size_t(1)}) {
        auto reject=[](auto&&...) -> T {assert(false);return 0;};
        auto full=lazy_multiply<T>(T(2),T(3),[](auto&&...)->std::pair<T,T>{assert(false);return {0,0};},n);
        big_vector<T>a(n,T(2));auto semi=lazy_multiply<T>(a,T(3),reject,n);
        auto e=lazy_exp<T>(T(0),reject,n);
        assert(full==big_vector<T>(n,T(6)) && semi==full && e==big_vector<T>(n,T(1)));
    }

}
