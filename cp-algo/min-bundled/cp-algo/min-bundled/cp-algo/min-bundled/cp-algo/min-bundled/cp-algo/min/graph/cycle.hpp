#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/min/graph/cycle.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/min/graph/cycle.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/min/graph/cycle.hpp"
#line 1 "cp-algo/min/graph/cycle.hpp"
#line 1 "cp-algo/min/graph/dfs.hpp"
#line 1 "cp-algo/min/graph/base.hpp"
#line 1 "cp-algo/min/graph/edge_types.hpp"
#include <iostream>
#include <cstdint>
namespace cp_algo::graph{using node_index=int;
struct edge_base{int xor_nodes;
edge_base(){}
edge_base(node_index from,node_index to):xor_nodes(from^to){}
node_index traverse(node_index from)const{return xor_nodes^from;}
static auto read(node_index v0=0){node_index u,v;
std::cin>>u>>v;
u-=v0;
v-=v0;
return std::pair{u,edge_base(u,v)};}};
struct weighted_edge:edge_base{int64_t w;
weighted_edge(){}
weighted_edge(node_index from,node_index to,int64_t w):edge_base(from,to),w(w){}
static auto read(node_index v0=0){auto[u,e]=edge_base::read(v0);
int64_t w;
std::cin>>w;
return std::pair{u,weighted_edge(u,e.traverse(u),w)};}};
template<typename edge>
concept edge_type=std::is_base_of_v<edge_base,edge>;
template<typename edge>
concept weighted_edge_type=std::is_base_of_v<weighted_edge,edge>;}
#line 1 "cp-algo/min/graph/concepts.hpp"
#line 4 "cp-algo/min/graph/concepts.hpp"
#include <type_traits>
namespace cp_algo::graph{enum graph_mode{directed,undirected};
template<typename T,typename=void>
struct graph_traits:std::false_type{};
template<typename T>
struct graph_traits<T,std::void_t<typename T::edge_t,decltype(T::mode)>>:std::true_type{using edge_t=typename T::edge_t;
static constexpr auto mode=T::mode;
static constexpr bool is_directed=mode==directed;
static constexpr bool is_undirected=mode==undirected;
static constexpr bool is_weighted=weighted_edge_type<edge_t>;};
template<typename G>
concept graph_type=graph_traits<G>::value;
template<typename G>
concept digraph_type=graph_type<G>&&graph_traits<G>::is_directed;
template<typename G>
concept undirected_graph_type=graph_type<G>&&graph_traits<G>::is_undirected;
template<typename G>
concept weighted_graph_type=graph_type<G>&&graph_traits<G>::is_weighted;
template<typename G>
concept weighted_digraph_type=digraph_type<G>&&graph_traits<G>::is_weighted;
template<typename G>
concept weighted_undirected_graph_type=undirected_graph_type<G>&&graph_traits<G>::is_weighted;}
#line 1 "cp-algo/min/structures/stack_union.hpp"
#line 1 "cp-algo/min/util/big_alloc.hpp"
#include <vector>
#include <cstddef>
#line 6 "cp-algo/min/util/big_alloc.hpp"
#if defined(__linux__) || defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#  define CP_ALGO_USE_MMAP 1
#  include <sys/mman.h>
#else
#  define CP_ALGO_USE_MMAP 0
namespace cp_algo{template<typename T,std::size_t Align=32>
class big_alloc{static_assert(Align>=alignof(void*),"Align must be at least pointer-size");
static_assert(std::popcount(Align)==1,"Align must be a power of two");
public:using value_type=T;
template<class U>struct rebind{using other=big_alloc<U,Align>;};
constexpr bool operator==(const big_alloc&)const=default;
constexpr bool operator!=(const big_alloc&)const=default;
big_alloc()noexcept=default;
template<typename U,std::size_t A>
big_alloc(const big_alloc<U,A>&)noexcept{}
[[nodiscard]]T*allocate(std::size_t n){std::size_t padded=round_up(n*sizeof(T));
std::size_t align=std::max<std::size_t>(alignof(T),Align);
#if CP_ALGO_USE_MMAP
if(padded>=MEGABYTE){void*raw=mmap(nullptr,padded,
PROT_READ|PROT_WRITE,
MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
madvise(raw,padded,MADV_HUGEPAGE);
madvise(raw,padded,MADV_POPULATE_WRITE);
return static_cast<T*>(raw);}
return static_cast<T*>(::operator new(padded,std::align_val_t(align)));}
void deallocate(T*p,std::size_t n)noexcept{if(!p)return;
std::size_t padded=round_up(n*sizeof(T));
std::size_t align=std::max<std::size_t>(alignof(T),Align);
#if CP_ALGO_USE_MMAP
if(padded>=MEGABYTE){munmap(p,padded);return;}::operator delete(p,padded,std::align_val_t(align));}
private:static constexpr std::size_t MEGABYTE=1<<20;
static constexpr std::size_t round_up(std::size_t x)noexcept{return(x+Align-1)/Align*Align;}};
template<typename T>
using big_vector=std::vector<T,big_alloc<T>>;}
#endif
#endif
#endif
#line 5 "cp-algo/min/structures/stack_union.hpp"
#include <iterator>
#include <ranges>
namespace cp_algo::structures{template<class datatype>
struct stack_union{stack_union(int n=0):head(n),next(1),data(1){}
void push(int v,datatype const&vdata){next.push_back(head[v]);
head[v]=(int)std::size(next)-1;
data.push_back(vdata);}
template<typename... Args>
void emplace(int v,Args&&... vdata){next.push_back(head[v]);
head[v]=(int)std::size(next)-1;
data.emplace_back(std::forward<Args>(vdata)...);}
void reserve(int m){data.reserve(m);
next.reserve(m);}
size_t size()const{return std::size(head);}
size_t nodes()const{return std::size(data);}
template<typename Su>
struct _iterator{using value_type=std::conditional_t<std::is_const_v<Su>,const datatype,datatype>;
using difference_type=std::ptrdiff_t;
Su*su=nullptr;
int sv=0;
value_type&operator*()const{return su->data[sv];}
_iterator&operator++(){sv=su->next[sv];
return*this;}
_iterator operator++(int){auto tmp=*this;++*this;return tmp;}
friend bool operator==(_iterator const&it,std::default_sentinel_t){return it.sv==0;}};
using iterator=_iterator<stack_union<datatype>>;
using const_iterator=_iterator<const stack_union<datatype>>;
auto operator[](this auto&&self,int v){using Iter=_iterator<std::remove_reference_t<decltype(self)>>;
return std::ranges::subrange(Iter{&self,self.head[v]},std::default_sentinel);}
big_vector<int>head,next;
big_vector<datatype>data;};}
#line 7 "cp-algo/min/graph/base.hpp"
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
#line 4 "cp-algo/min/graph/dfs.hpp"
#include <variant>
#include <stack>
namespace cp_algo::graph{enum node_state{unvisited,visiting,visited,blocked};
template<graph_type graph>
struct dfs_context{big_vector<node_state>state;
graph const*g;
bool done=false;
dfs_context(graph const&g):state(g.n()),g(&g){}
void on_enter(node_index){}
void on_tree_edge(node_index,edge_index){}
void on_return_from_child(node_index,edge_index){}
void on_back_edge(node_index,edge_index){}
void on_forward_cross_edge(node_index,edge_index){}
void on_exit(node_index){}};
template<template<typename>class Context,graph_type graph>
Context<graph>dfs(graph const&g){Context<graph>context(g);
auto const&adj=g.incidence_lists();
struct frame{node_index v;
[[no_unique_address]]std::conditional_t<
undirected_graph_type<graph>,
edge_index,std::monostate>ep;
int sv;
enum{INIT,PROCESS_EDGES,HANDLE_CHILD}state;};
std::stack<frame>dfs_stack;
for(auto root:g.nodes()){if(context.done)break;
if(context.state[root]!=unvisited)continue;
if constexpr(undirected_graph_type<graph>){dfs_stack.push({root,-1,0,frame::INIT});}else{dfs_stack.push({root,{},0,frame::INIT});}
while(!empty(dfs_stack)){auto&f=dfs_stack.top();
if(f.state==frame::INIT){context.state[f.v]=visiting;
context.on_enter(f.v);
f.sv=adj.head[f.v];
f.state=frame::PROCESS_EDGES;
continue;}
if(f.state==frame::HANDLE_CHILD){auto e=adj.data[f.sv];
f.sv=adj.next[f.sv];
context.on_return_from_child(f.v,e);
f.state=frame::PROCESS_EDGES;
continue;}
bool found_child=false;
while(f.sv!=0&&!context.done){auto e=adj.data[f.sv];
if constexpr(undirected_graph_type<graph>){if(f.ep==e){f.sv=adj.next[f.sv];
continue;}}
node_index u=g.edge(e).traverse(f.v);
if(context.state[u]==unvisited){context.on_tree_edge(f.v,e);
f.state=frame::HANDLE_CHILD;
if constexpr(undirected_graph_type<graph>){dfs_stack.push({u,e,0,frame::INIT});}else{dfs_stack.push({u,{},0,frame::INIT});}
found_child=true;
break;}else if(context.state[u]==visiting){context.on_back_edge(f.v,e);}else if(context.state[u]==visited){context.on_forward_cross_edge(f.v,e);}
f.sv=adj.next[f.sv];}
if(found_child)continue;
context.state[f.v]=visited;
context.on_exit(f.v);
dfs_stack.pop();}}
return context;}}
#line 5 "cp-algo/min/graph/cycle.hpp"
#include <deque>
namespace cp_algo::graph{template<graph_type graph>
struct cycle_context:dfs_context<graph>{using base=dfs_context<graph>;
using base::base;
std::deque<edge_index>cycle;
bool closed=false;
int v0;
void on_return_from_child(node_index v,edge_index e){if(!empty(cycle)&&!closed){cycle.push_front(e);
closed|=v==v0;}}
void on_back_edge(node_index v,edge_index e){if(empty(cycle)){v0=base::g->edge(e).traverse(v);
base::done=true;
closed=v==v0;
cycle.push_front(e);}}};
template<graph_type graph>
std::pair<node_index,std::deque<edge_index>>find_cycle(graph const&g){auto context=dfs<cycle_context>(g);
return{context.v0,context.cycle};}}