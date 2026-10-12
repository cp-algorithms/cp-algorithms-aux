#include "cp-algo/linalg/matrix.hpp"
#include <cassert>
#include <random>
using namespace cp_algo::math; using namespace cp_algo::linalg;
template<class T> void check() {
    std::mt19937 rng(17);
    for(int n:{1,2,4,5,9,17}) {
        matrix<T>a(n,n),b(n,n);
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){a[i][j]=rng()%2?T(T::mod()-1):T(rng());b[i][j]=rng()%2?T(T::mod()-1):T(rng());}
        auto c=a*b;
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){T want=0;for(int k=0;k<n;k++)want+=a[i][k]*b[k][j];assert(c[i][j]==want);}
    }
    matrix<T>a(1,4),b(4,1);for(int i=0;i<4;i++)a[0][i]=b[i][0]=T(T::mod()-1);assert((a*b)[0][0]==T(4));
}
int main(){check<modint<2147483647LL>>();check<modint<2013265921LL>>();check<modint<2305843009213693951LL>>();dynamic_modint<int64_t>::with_mod(2147483647,[]{check<dynamic_modint<int64_t>>();});}
