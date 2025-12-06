#line 1 "cp-algo/min-bundled/cp-algo/structures/dsu.hpp"
#line 1 "cp-algo/structures/dsu.hpp"
#include <numeric>
#include <vector>
namespace cp_algo::structures{struct disjoint_set_union{disjoint_set_union(int n):par(n){std::iota(begin(par),end(par),0);}
int get(int v){return v==par[v]?v:par[v]=get(par[v]);}
bool uni(int a,int b){a=get(a);
b=get(b);
par[a]=b;
return a!=b;}
private:std::vector<int>par;};
using dsu=disjoint_set_union;}