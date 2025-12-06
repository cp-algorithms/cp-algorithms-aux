#ifndef CP_ALGO_GRAPH_BASE_HPP
#define CP_ALGO_GRAPH_BASE_HPP
#include "edge_types.hpp"
#include "concepts.hpp"
#include "../structures/stack_union.hpp"
#include <ranges>
namespace cp_algo::graph{using edge_index=int;
template<edge_type _edge_t=edge_base,graph_mode _mode=undirected>
struct graph{using edge_t=_edge_t;
static constexpr auto mode=_mode;
using incidence_list=structures::stack_union<edge_index>;
graph(int n,int v0=0):v0(v0),adj(n){}
graph transpose()const{static_assert(mode==directed,"transpose is only defined for directed graphs");
graph<edge_t,mode>gt(n(),v0);
for(auto v:nodes()){for(auto e:outgoing(v)){gt.add_edge(edge(e).traverse(v),edge(e));}}
return gt;}
edge_index add_edge(node_index u,edge_t e){edge_index idx=(edge_index)size(E);
E.push_back(e);
adj.push(u,idx);
if constexpr(mode==undirected){adj.push(e.traverse(u),idx);}
return idx;}
edge_index add_edge(node_index u,auto... Args){return add_edge(u,edge_t(u,Args...));}
void read_edges(node_index m){adj.reserve(mode==undirected?2*m:m);
for(edge_index i=0;i<m;i++){auto[u,e]=edge_t::read(v0);
add_edge(u,e);}}
auto outgoing(node_index v)const{return adj[v];}
auto edges()const{return E|std::views::all;}
auto nodes()const{return std::views::iota(node_index(0),n());}
auto edge_indices()const{return std::views::iota(edge_index(0),m());}
auto&&incidence_lists(this auto&&self){return self.adj;}
auto&&edge(this auto&&self,edge_index e){return self.E[e];}
node_index n()const{return(node_index)incidence_lists().size();}
edge_index m()const{return(edge_index)edges().size();}
private:node_index v0;
big_vector<edge_t>E;
incidence_list adj;};
template<edge_type edge_t=edge_base>
using digraph=graph<edge_t,directed>;
template<weighted_edge_type edge_t=weighted_edge,graph_mode mode=undirected>
using weighted_graph=graph<edge_t,mode>;
template<weighted_edge_type edge_t=weighted_edge>
using weighted_digraph=digraph<edge_t>;}
#endif