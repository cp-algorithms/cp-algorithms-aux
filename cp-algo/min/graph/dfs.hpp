#ifndef CP_ALGO_GRAPH_DFS_HPP
#define CP_ALGO_GRAPH_DFS_HPP
#include "base.hpp"
#include <variant>
#include <stack>
namespace cp_algo::graph{
enum node_state{unvisited,visiting,visited,blocked};
template<graph_type graph>
struct dfs_context{
big_vector<node_state>state;
graph const*g;
bool done=false;
dfs_context(graph const&g):state(g.n()),g(&g){}
void on_enter(node_index){}
void on_tree_edge(node_index,edge_index){}
void on_return_from_child(node_index,edge_index){}
void on_back_edge(node_index,edge_index){}
void on_forward_cross_edge(node_index,edge_index){}
void on_exit(node_index){}
};
template<template<typename>class Context,graph_type graph>
Context<graph>dfs(graph const&g){
Context<graph>context(g);
auto const&adj=g.incidence_lists();
struct frame{
node_index v;
[[no_unique_address]]std::conditional_t<
undirected_graph_type<graph>,
edge_index,std::monostate>ep;
int sv;
enum{INIT,PROCESS_EDGES,HANDLE_CHILD}state;
};
std::stack<frame>dfs_stack;
for(auto root:g.nodes()){
if(context.done)break;
if(context.state[root]!=unvisited)continue;
if constexpr(undirected_graph_type<graph>){
dfs_stack.push({root,-1,0,frame::INIT});
}else{
dfs_stack.push({root,{},0,frame::INIT});
}
while(!empty(dfs_stack)){
auto&f=dfs_stack.top();
if(f.state==frame::INIT){
context.state[f.v]=visiting;
context.on_enter(f.v);
f.sv=adj.head[f.v];
f.state=frame::PROCESS_EDGES;
continue;
}
if(f.state==frame::HANDLE_CHILD){
auto e=adj.data[f.sv];
f.sv=adj.next[f.sv];
context.on_return_from_child(f.v,e);
f.state=frame::PROCESS_EDGES;
continue;
}
bool found_child=false;
while(f.sv!=0&&!context.done){
auto e=adj.data[f.sv];
if constexpr(undirected_graph_type<graph>){
if(f.ep==e){
f.sv=adj.next[f.sv];
continue;
}
}
node_index u=g.edge(e).traverse(f.v);
if(context.state[u]==unvisited){
context.on_tree_edge(f.v,e);
f.state=frame::HANDLE_CHILD;
if constexpr(undirected_graph_type<graph>){
dfs_stack.push({u,e,0,frame::INIT});
}else{
dfs_stack.push({u,{},0,frame::INIT});
}
found_child=true;
break;
}else if(context.state[u]==visiting){
context.on_back_edge(f.v,e);
}else if(context.state[u]==visited){
context.on_forward_cross_edge(f.v,e);
}
f.sv=adj.next[f.sv];
}
if(found_child)continue;
context.state[f.v]=visited;
context.on_exit(f.v);
dfs_stack.pop();
}
}
return context;
}
}
#endif
