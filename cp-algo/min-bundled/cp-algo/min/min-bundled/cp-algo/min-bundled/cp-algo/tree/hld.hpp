#line 1 "cp-algo/min/min-bundled/cp-algo/min-bundled/cp-algo/tree/hld.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/tree/hld.hpp"
#line 1 "cp-algo/tree/hld.hpp"
#line 1 "cp-algo/tree/ascending_dfs.hpp"
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
#line 5 "cp-algo/tree/ascending_dfs.hpp"
#include <cassert>
#line 8 "cp-algo/tree/ascending_dfs.hpp"
namespace cp_algo::graph{template<undirected_graph_type graph>
void ascending_dfs(graph const&tree,auto&degree,auto&&next,auto&&callback,node_index root){for(auto v:tree.nodes()){while(degree[v]==1){edge_index ep=next(v);
callback(v,ep);
degree[v]--;
v=tree.edge(ep).traverse(v);
degree[v]--;}}
callback(root,-1);}
template<undirected_graph_type graph>
auto xor_dfs(graph const&tree,auto&&callback,node_index root=0){std::vector<edge_index>neig_xor(tree.n());
std::vector<int>degree(tree.n());
for(auto v:tree.nodes()){degree[v]=(int)std::ranges::distance(tree.outgoing(v));}
degree[root]=0;
for(auto v:tree.nodes()){for(auto e:tree.outgoing(v)){neig_xor[v]^=e;}}
neig_xor[root]^=edge_index(-1);
ascending_dfs(tree,degree,[&](auto v){edge_index ep=neig_xor[v];
neig_xor[tree.edge(ep).traverse(v)]^=ep;
return ep;},callback,root);
return neig_xor;}
template<undirected_graph_type graph>
void parent_dfs(graph const&tree,std::vector<edge_index>const&parent,auto&&callback){std::vector<int>degree(tree.n());
node_index root=-1;
for(auto[v,e]:parent|std::views::enumerate){if(e!=-1){degree[v]++;
degree[tree.edge(e).traverse(node_index(v))]++;}else{root=node_index(v);}}
assert(root!=-1);
degree[root]=0;
ascending_dfs(tree,degree,[&](auto v){return parent[v];},callback,root);}}
#line 6 "cp-algo/tree/hld.hpp"
namespace cp_algo::graph{struct heavy_light{big_vector<node_index>size,in,up,par;
template<undirected_graph_type graph>
heavy_light(graph const&g,node_index root=0,std::vector<edge_index>const*parents_ptr=nullptr):size(g.n(),1),in(g.n()),up(g.n()),par(g.n()){big_vector<node_index>topsort;
topsort.reserve(g.n());
auto push_size=[&](node_index v,edge_index e){if(size[v]>1){topsort.push_back(v);}
if(v!=root){auto p=g.edge(e).traverse(v);
size[p]+=size[v];
par[v]=p;}};
if(parents_ptr){parent_dfs(g,*parents_ptr,push_size);}else{xor_dfs(g,push_size,root);}
par[root]=up[root]=root;
for(auto v:topsort|std::views::reverse){node_index big=-1;
for(auto e:g.outgoing(v)){auto u=g.edge(e).traverse(v);
if(size[u]>size[v])continue;
if(big==-1||size[u]>size[big]){big=u;}}
int t=in[v]+size[big];
for(auto e:g.outgoing(v)){auto u=g.edge(e).traverse(v);
if(size[u]>size[v])continue;
if(u==big){in[u]=in[v]+1;
up[u]=up[v];}else{in[u]=t+1;
t+=size[u];
up[u]=u;}}}}
enum lca_mode{without_distances,with_distances};
template<lca_mode mode=without_distances>
auto lca(node_index a,node_index b){int dista=0,distb=0;
while(up[a]!=up[b]){if(in[up[a]]<in[up[b]]){if constexpr(mode==with_distances)distb+=in[b]-in[up[b]]+1;
b=par[up[b]];}else{if constexpr(mode==with_distances)dista+=in[a]-in[up[a]]+1;
a=par[up[a]];}}
node_index c=in[a]<in[b]?a:b;
if constexpr(mode==with_distances){return std::tuple{c,dista+in[a]-in[c],distb+in[b]-in[c]};}else{return c;}}
big_vector<node_index>rin;
void compute_rin(){if(empty(rin)){rin.resize(std::size(in));
for(auto[v,inv]:in|std::views::enumerate){rin[inv]=node_index(v);}}}
node_index jump_up(node_index v,int steps){compute_rin();
while(steps>0){int path_dist=in[v]-in[up[v]];
if(steps<=path_dist){return rin[in[v]-steps];}
steps-=path_dist+1;
v=par[up[v]];}
return v;}
std::optional<node_index>jump(node_index from,node_index to,int steps){compute_rin();
auto[l,dist_from,dist_to]=lca<with_distances>(from,to);
auto dist=dist_from+dist_to;
if(steps>dist)return std::nullopt;
if(steps<=dist_from){return jump_up(from,steps);}else{return jump_up(to,dist-steps);}}};}