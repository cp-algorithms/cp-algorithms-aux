#ifndef CP_ALGO_STRUCTURES_DSU_HPP
#define CP_ALGO_STRUCTURES_DSU_HPP
#include "../util/big_alloc.hpp"
#include <numeric>
#include <vector>
namespace cp_algo::structures{struct disjoint_set_union{disjoint_set_union(int n):par(n){std::iota(begin(par),end(par),0);}int get(int v){return v==par[v]?v:par[v]=get(par[v]);}bool uni(int a,int b){a=get(a);b=get(b);par[a]=b;return a!=b;}private:big_vector<int>par;};using dsu=disjoint_set_union;}
#endif
