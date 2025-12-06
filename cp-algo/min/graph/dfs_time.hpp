#ifndef CP_ALGO_GRAPH_DFS_TIME_HPP
#define CP_ALGO_GRAPH_DFS_TIME_HPP
#include "dfs.hpp"
#include "base.hpp"
namespace cp_algo::graph{
template<graph_type graph>
struct dfs_time_context:dfs_context<graph>{
using base=dfs_context<graph>;
big_vector<int>tin;
int timer;
dfs_time_context(graph const&g):base(g),tin(g.n()),timer(0){}
void on_enter(node_index v){
tin[v]=timer++;
}
};
template<graph_type graph>
struct dfs_time_range_context:dfs_time_context<graph>{
using base=dfs_time_context<graph>;
big_vector<int>tout;
dfs_time_range_context(graph const&g):base(g),tout(g.n()){}
void on_exit(node_index v){
tout[v]=base::timer;
}
};
template<graph_type graph>
struct dfs_low_context:dfs_time_context<graph>{
using base=dfs_time_context<graph>;
big_vector<int>low;
dfs_low_context(graph const&g):base(g),low(g.n()){}
void on_enter(node_index v){
base::on_enter(v);
low[v]=base::tin[v];
}
void on_return_from_child(node_index v,edge_index e){
node_index u=base::g->edge(e).traverse(v);
low[v]=std::min(low[v],low[u]);
}
void on_back_edge(node_index v,edge_index e){
node_index u=base::g->edge(e).traverse(v);
low[v]=std::min(low[v],base::tin[u]);
}
void on_forward_cross_edge(node_index v,edge_index e){
node_index u=base::g->edge(e).traverse(v);
low[v]=std::min(low[v],base::tin[u]);
}
};
}
#endif
