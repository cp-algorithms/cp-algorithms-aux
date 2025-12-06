#ifndef CP_ALGO_GRAPH_MST_HPP
#define CP_ALGO_GRAPH_MST_HPP
#include "base.hpp"
#include "../structures/dsu.hpp"
#include "../util/sort.hpp"
#include <algorithm>
namespace cp_algo::graph{
template<weighted_undirected_graph_type graph>
std::pair<int64_t,std::vector<edge_index>>mst(graph const&g){
struct edge{
edge_index idx;
node_index v;
};
std::vector<edge>edges;
for(auto v:g.nodes()){
for(auto e:g.outgoing(v)){
if(v<g.edge(e).traverse(v)){
edges.push_back({e,v});
}
}
}
radix_sort(edges,[&](auto e){
return g.edge(e.idx).w;
});
structures::dsu me(g.n());
int64_t total=0;
std::vector<edge_index>mst;
for(auto[idx,v]:edges){
if(me.uni(v,g.edge(idx).traverse(v))){
total+=g.edge(idx).w;
mst.push_back(idx);
}
}
return{total,mst};
}
}
#endif
