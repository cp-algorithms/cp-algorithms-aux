#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/graph/shortest_path.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/graph/shortest_path.hpp"
#line 1 "cp-algo/graph/shortest_path.hpp"
#line 1 "cp-algo/graph/base.hpp"
#line 1 "cp-algo/graph/edge_types.hpp"
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
#line 1 "cp-algo/graph/concepts.hpp"
#line 4 "cp-algo/graph/concepts.hpp"
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
#line 1 "cp-algo/structures/stack_union.hpp"
#line 1 "cp-algo/util/big_alloc.hpp"
#include <vector>
#include <cstddef>
#line 7 "cp-algo/util/big_alloc.hpp"
#if defined(__linux__) || defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#  define CP_ALGO_USE_MMAP 1
#  include <sys/mman.h>
#else
#  define CP_ALGO_USE_MMAP 0
#endif
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
#endif
return static_cast<T*>(::operator new(padded,std::align_val_t(align)));}
void deallocate(T*p,std::size_t n)noexcept{if(!p)return;
std::size_t padded=round_up(n*sizeof(T));
std::size_t align=std::max<std::size_t>(alignof(T),Align);
#if CP_ALGO_USE_MMAP
if(padded>=MEGABYTE){munmap(p,padded);return;}
#endif::operator delete(p,padded,std::align_val_t(align));}
private:static constexpr std::size_t MEGABYTE=1<<20;
static constexpr std::size_t round_up(std::size_t x)noexcept{return(x+Align-1)/Align*Align;}};
template<typename T>
using big_vector=std::vector<T,big_alloc<T>>;}
#line 5 "cp-algo/structures/stack_union.hpp"
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
#line 7 "cp-algo/graph/base.hpp"
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
#line 4 "cp-algo/graph/shortest_path.hpp"
#include <algorithm>
#include <queue>
namespace cp_algo::graph{struct shortest_path_context{std::vector<int64_t>dist;
std::vector<edge_index>pre;
static constexpr int64_t inf=1e18;
shortest_path_context(int n):dist(n,inf),pre(n){}};
struct dijkstra_context:shortest_path_context{struct que_t{int64_t dist;
node_index v;
bool operator<(que_t const&other)const{return dist>other.dist;}};
std::priority_queue<que_t>pq;
dijkstra_context(int n):shortest_path_context(n){}
void push(node_index,edge_index,node_index v){pq.push({dist[v],v});}
std::optional<node_index>next_node(){while(!empty(pq)){auto[dv,v]=pq.top();
pq.pop();
if(dv==dist[v]){return v;}}
return std::nullopt;}};
struct spfa_context:shortest_path_context{std::queue<node_index>que;
std::vector<char>flags;
static constexpr char in_queue=1;
static constexpr char invalidated=2;
spfa_context(int n):shortest_path_context(n),flags(n){}
void push(node_index,edge_index,node_index v){if(!(flags[v]&in_queue)){que.push(v);
flags[v]|=in_queue;}}
std::optional<node_index>next_node(){while(!que.empty()){node_index v=que.front();
que.pop();
flags[v]&=~in_queue;
if(!(flags[v]&invalidated)){return v;}}
return std::nullopt;}};
struct deep_spfa_context:spfa_context{struct traverse_edge{edge_index e;
node_index v;};
std::vector<std::basic_string<traverse_edge>>dependents;
deep_spfa_context(int n):spfa_context(n),dependents(n){}
void push(node_index u,edge_index e,node_index v){invalidate_subtree(v);
dependents[u].push_back({e,v});
flags[v]&=~invalidated;
spfa_context::push(u,e,v);}
void invalidate_subtree(node_index v){std::vector<node_index>to_invalidate={v};
while(!empty(to_invalidate)){node_index u=to_invalidate.back();
to_invalidate.pop_back();
flags[u]|=invalidated;
flags[u]&=~in_queue;
for(auto[e,v]:dependents[u]){if(pre[v]==e){to_invalidate.push_back(v);}}
dependents[u].clear();}}};
template<typename Context,weighted_graph_type graph>
Context sssp_impl(graph const&g,node_index s){Context context(g.n());
context.dist[s]=0;
context.pre[s]=-1;
context.push(s,-1,s);
while(auto ov=context.next_node()){node_index v=*ov;
for(auto e:g.outgoing(v)){node_index u=g.edge(e).traverse(v);
auto w=g.edge(e).w;
if(context.dist[v]+w<context.dist[u]){context.dist[u]=context.dist[v]+w;
context.pre[u]=e;
context.push(v,e,u);}}}
return context;}
template<weighted_graph_type graph>
shortest_path_context dijkstra(graph const&g,node_index s){return sssp_impl<dijkstra_context>(g,s);}
template<weighted_graph_type graph>
shortest_path_context spfa(graph const&g,node_index s){return sssp_impl<spfa_context>(g,s);}
template<weighted_graph_type graph>
shortest_path_context deep_spfa(graph const&g,node_index s){return sssp_impl<deep_spfa_context>(g,s);}
template<weighted_graph_type graph>
shortest_path_context single_source_shortest_path(graph const&g,node_index s){bool negative_edges=false;
for(auto e:g.edges()){negative_edges|=e.w<0;}
return negative_edges?deep_spfa(g,s):dijkstra(g,s);}
std::vector<edge_index>recover_path(auto const&g,auto const&pre,node_index s,node_index t){std::vector<edge_index>path;
node_index v=t;
while(v!=s){path.push_back(pre[v]);
v=g.edge(pre[v]).traverse(v);}
std::ranges::reverse(path);
return path;}
template<weighted_graph_type graph>
std::optional<std::pair<int64_t,std::vector<edge_index>>>shortest_path(graph const&g,node_index s,node_index t){auto[dist,pre]=single_source_shortest_path(g,s);
if(dist[t]==shortest_path_context::inf){return std::nullopt;}
return{{dist[t],recover_path(g,pre,s,t)}};}}