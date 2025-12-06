#line 1 "verify/tree/jump.test.cpp"
#define PROBLEM "https:#pragma GCC optimize("Ofast,unroll-loops")
#include <iostream>
#line 1 "blazingio/blazingio.min.hpp"
#define M$(x,...)_mm256_##x##_epi8(__VA_ARGS__)
#define $u(...)__VA_ARGS__
#if __APPLE__
#define $m(A,B)A
#else
#define $m(A,B)B
#endif
#if _WIN32
#define $w(A,B)A
#else
#define $w(A,B)B
#endif
#if __i386__|_M_IX86
#define $H(A,B)A
#else
#define $H(A,B)B
#endif
#if __aarch64__
#define $a(A,B)A
#else
#define $a(A,B)B
#endif
#define $P(x)void F(x K){#define $T template<$c T
#define $c class
#define $C constexpr
#define $R return
#define $O operator
#define u$ uint64_t
#define $r $R*this;
#include<array>
#include<bitset>
#include<complex>
#include<cstring>
#include $a(<arm_neon.h>,<immintrin.h>)
#include<stdint.h>
#include $w(<windows.h>,<sys/mman.h>)
#include<sys/stat.h>
#include $w(<io.h>,<unistd.h>)
#include $w(<ios>,<sys/resource.h>)
#if _MSC_VER
#define __builtin_add_overflow(a,b,c)_addcarry_u64(0,a,b,c)
#define $s
#else
$H(,u$ _umul128(u$ a,u$ b,u$*D){auto x=(__uint128_t)a*b;*D=u$(x>>64);$R(u$)x;})
#define $s $a(,__attribute__((target("avx2"))))
#endif
#define $z $a(16,32)
#define $t $a(uint8x16_t,__m256i)
#define $I $w(__forceinline,__attribute__((always_inline)))
#define $F M(),
#define E$(x)if(!(x))abort();
$w(LONG WINAPI $x(_EXCEPTION_POINTERS*);,)namespace $f{using namespace std;struct B{enum $c A:char{}c;B&$O=(char x){c=A{x};$r}$O char(){$R(char)c;}};$C u$ C=~0ULL/255;struct D{string&K;};static B E[65568];template<int F>struct G{B*H,*S;void K(off_t C){$w(char*D=(char*)VirtualAlloc(0,(C+8191)&-4096,8192,1);E$(D)E$(VirtualFree(D,0,32768))DWORD A=C&-65536;E$(!A||MapViewOfFileEx(CreateFileMapping(GetStdHandle(-10),0,2,0,A,0),4,0,0,0,D)==D)E$(VirtualAlloc(D+A,65536,12288,4)==D+A)E$(~_lseek(0,A,0))DWORD E=0;ReadFile(GetStdHandle(-10),D+A,65536,&E,0);,int A=getpagesize();char*D=(char*)mmap(0,C+A,3,2,0,0);E$(D!=(void*)-1)E$(mmap(D+((C+A-1)&-A),A,3,$m(4114,50),-1,0)!=(void*)-1))H=(B*)D+C;*H=10;H[1]=48;H[2]=0;S=(B*)D;}void L(){H=S=E;}$I void M(){if(F&&S==H){$w(DWORD A=0;ReadFile(GetStdHandle(-10),S=E,65536,&A,0);,$a($u(register long A asm("x0")=0,D asm("x1")=(long)E,G asm("x2")=65536,C asm($m("x16","x8"))=$m(3,63);asm volatile("svc 0"$m("x80",):"+r"(A),"+r"(D):"r"(C),"r"(G));S=launder(E);),off_t A=$H(3,$m(33554435,0));B*D=E;asm volatile($H("int $128","syscall"):"+a"(A),$H("+c"(D):"b","+S"(D):"D")(0),"d"(65536)$H(,$u(:"rcx","r11")));S=D;))H=S+A;*H=10;if(!A)E[1]=48,E[2]=0;}}$T>$I void N(T&x){while($F(*S&240)==48)x=T(x*10+(*S++-48));}$T>$I decltype((void)~T{1})O(T&x){M();int A=is_signed_v<T>&&*S==45;S+=A;N(x=0);x=A?1+~x:x;}$T>$I decltype((void)T{1.})O(T&x){M();int A=*S==45;S+=A;$F S+=*S==43;u$ n=0;int i=0;for(;i<18&&($F*S&240)==48;i++)n=n*10+*S++-48;int B=20;int C=*S==46;S+=C;for(;i<18&&($F*S&240)==48;i++)n=n*10+*S++-48,B-=C;x=(T)n;while(($F*S&240)==48)x=x*10+*S++-48,B-=C;if(*S==46)S++,C=1;while(($F*S&240)==48)x=x*10+*S++-48,B-=C;int D;if((*S|32)==101)S++,$F S+=*S==43,O(D),B+=D;static $C auto E=[](){array<T,41>E{};T x=1;for(int i=21;i--;)E[40-i]=x,E[i]=1/x,x*=10;$R E;}();while(B>40)x*=(T)1e10,B-=10;while(B<0)x*=(T)1e-10,B+=10;x*=E[B];x=A?-x:x;}$I void O(bool&x){$F x=*S++==49;}$I void O(char&x){$F x=*S++;}$I void O(uint8_t&x){$F x=*S++;}$I void O(int8_t&x){$F x=*S++;}$T>$s void P(string&K,T C){M();B*G=S;C();K.assign((char*)G,S-G);while(F&&S==H&&($F H!=E)){C();K.append(E,S);}}$s void O(string&K){P(K,[&]()$s{B*p=S;$w(ULONG R;,)$t x;$a(uint64x2_t A;while(memcpy(&x,p,16),A=uint64x2_t(x<33),!(A[0]|A[1]))p+=16;S=p+(A[0]?0:8)+$w((_BitScanForward64(&R,A[0]?A[0]:A[1]),R),__builtin_ctzll(A[0]?A[0]:A[1]))/8;,int J;$t C=M$(set1,32);while(memcpy(&x,p,32),!(J=M$(movemask,M$(cmpeq,C,_mm256_max_epu8(C,x)))))p+=32;S=p+$w((_BitScanForward(&R,J),R),__builtin_ctz(J));)});}$s void O(D&A){P(A.K,[&](){S=(B*)memchr(S,10,H-S+1);});if(A.K.size()&&A.K.back()==13)A.K.pop_back();if(A.K.empty()||S<H)S+=*S==10;}$T>$I void O(complex<T>&K){T A,B{};if($F*S==40){S++;O(A);if($F*S++==44)Q(B),S++;}else O(A);K={A,B};}template<size_t N>$s void O(bitset<N>&K){if(N>4095&&!*this)$R;ptrdiff_t i=N;while(i)if($F i%$z||H-S<$z)K[--i]=*S++==49;else{B*p=S;for(int64_t j=0;j<min(i,H-S)/$z;j++){i-=$z;$t x;memcpy(&x,p,$z);$a(auto B=(uint8x16_t)vdupq_n_u64(~2ULL/254)&(48-x);auto C=vzip_u8(vget_high_u8(B),vget_low_u8(B));auto y=vaddvq_u16((uint16x8_t)vcombine_u8(C.val[0],C.val[1]));,u$ a=~0ULL/65025;auto y=$w(_byteswap_ulong,__builtin_bswap32)(M$(movemask,M$(shuffle,_mm256_slli_epi32(x,7),_mm256_set_epi64x(a+C*24,a+C*16,a+C*8,a))));)p+=$z;memcpy((char*)&K+i/8,&y,$z/8);}S=p;}}$T>$I void Q(T&K){if(!is_same_v<T,D>)while($F(uint8_t)*S<33)S++;O(K);}$O bool(){$R!!*this;}bool $O!(){$R S>H;}};struct U{G<0>A;G<1>B;U(){struct stat D;E$(~fstat(0,&D))(D.st_mode>>12)==8?A.K(D.st_size):B.L();}U*tie(nullptr_t){$R this;}void sync_with_stdio(bool){}$T>$I U&$O>>(T&K){A.S?A.Q(K):B.Q(K);$r}$O bool(){$R!!*this;}bool $O!(){$R A.S?!A:!B;}};short A[100];char L[64]{1};struct
V{char*D;B*S;int J;V(){$w(E$(D=(char*)VirtualAlloc(0,536870912,8192,4))E$(VirtualAlloc(D,4096,4096,260))AddVectoredExceptionHandler(1,$x);,size_t C=536870912;$m(,rlimit E;getrlimit(RLIMIT_AS,&E);if(~E.rlim_cur)C=25165824;)D=(char*)mmap(0,C,3,$m(4162,16418),-1,0);E$(D!=(void*)-1))S=(B*)D;for(int i=0;i<100;i++)A[i]=short((48+i/10)|((48+i%10)<<8));for(int i=1;i<64;i++)L[i]=L[i-1]+(0x8922489224892249>>i&1);}~V(){flush($w(!J,));}void flush($w(int F=0,)){$w(J=1;auto E=GetStdHandle(-11);auto C=F?ReOpenFile(E,1073741824,7,2684354560):(void*)-1;DWORD A;E$(C==(void*)-1?WriteFile(E,D,DWORD((char*)S-D),&A,0):(WriteFile(C,D,DWORD(((char*)S-D+4095)&-4096),&A,0)&&~_chsize(1,int((char*)S-D)))),auto G=D;ssize_t A;while((A=write(1,G,(char*)S-G))>0)G+=A;E$(~A))S=(B*)D;}$P(char)*S++=K;}$P(uint8_t)*S++=K;}$P(int8_t)*S++=K;}$P(bool)*S++=48+K;}$T>decltype((void)~T{1})F(T K){using D=make_unsigned_t<T>;D C=K;if(K<0)F('-'),C=1+~C;static $C auto N=[](){array<D,5*sizeof(T)/2>N{};D n=1;for(size_t i=1;i<N.size();i++)n*=10,N[i]=n;$R N;}();$w(ULONG M;,)int G=L[$w(($H(_BitScanReverse(&M,ULONG((int64_t)C>>32))?M+=32:_BitScanReverse(&M,(ULONG)C|1),_BitScanReverse64(&M,C|1)),M),63^__builtin_clzll(C|1))];G-=C<N[G-1];short H[20];if $C(sizeof(T)==2){auto n=33555U*C-C/2;u$ H=A[n>>25];n=(n&33554431)*25;H|=A[n>>23]<<16;H|=u$(48+((n&8388607)*5>>22))<<32;H>>=40-G*8;memcpy(S,&H,8);}else if $C(sizeof(T)==4){auto n=1441151881ULL*C;$H(n>>=25;n++;for(int i=0;i<5;i++){H[i]=A[n>>32];n=(n&~0U)*100;},int K=57;auto J=~0ULL>>7;for(int i=0;i<5;i++){H[i]=A[n>>K];n=(n&J)*25;K-=2;J/=4;})memcpy(S,(B*)H+10-G,16);}else{$H($u(if(C<(1ULL<<32)){$R F((uint32_t)C);}auto J=(u$)1e10;auto x=C/J,y=C%J;int K=100000,b[]{int(x/K),int(x%K),int(y/K),int(y%K)};B H[40];for(int i=0;i<4;i++){int n=int((429497ULL*b[i]>>7)+1);B*p=H+i*5;*p=48+char(n>>25);n=(n&~0U>>7)*25;memcpy(p+1,A+(n>>23),2);memcpy(p+3,A+((n&~0U>>9)*25>>21),2);}),$u(u$ D,E=_umul128(18,C,&D),F;_umul128(0x725dd1d243aba0e8,C,&F);D+=__builtin_add_overflow(E,F+1,&E);for(int i=0;i<10;i++)H[i]=A[D],E=_umul128(100,E,&D);))memcpy(S,(B*)H+20-G,20);}S+=G;}$T>decltype((void)T{1.})F(T K){if(K<0)F('-'),K=-K;auto G=[&](){auto x=u$(K*1e12);$H($u(x-=x>999999999999;uint32_t n[]{uint32_t(x/1000000*429497>>7)+1,uint32_t(x%1000000*429497>>7)+1};int K=25,J=~0U>>7;for(int i=0;i<3;i++){for(int j=0;j<2;j++)memcpy(S+i*2+j*6,A+(n[j]>>K),2),n[j]=(n[j]&J)*25;K-=2;J/=4;}S+=12;),$u(u$ D,E=_umul128(472236648287,x,&D)>>8;E|=D<<56;D>>=8;E++;for(int i=0;i<6;i++)memcpy(S,A+D,2),S+=2,E=_umul128(100,E,&D);))};if(K==0)$R F('0');if(K>=1e16){K*=(T)1e-16;int B=16;while(K>=1)K*=(T).1,B++;F("0.");G();F('e');F(B);}else if(K>=1){auto B=(u$)K;F(B);if((K-=(T)B)>0)F('.'),G();}else F("0."),G();}$P(const char*)$w(size_t A=strlen(K);memcpy((char*)S,K,A);S+=A;,S=(B*)stpcpy((char*)S,K);)}$P(const uint8_t*)F((char*)K);}$P(const int8_t*)F((char*)K);}$P(string_view)memcpy(S,K.data(),K.size());S+=K.size();}$T>$P(complex<T>)*this<<'('<<K.real()<<','<<K.imag()<<')';}template<size_t N>$s $P(const bitset<N>&)auto i=N;while(i%$z)*S++=48+K[--i];B*p=S;while(i){i-=$z;$a(short,int)x;memcpy(&x,(char*)&K+i/8,$z/8);$a(auto A=(uint8x8_t)vdup_n_u16(x);vst1q_u8((uint8_t*)p,48-vtstq_u8(vcombine_u8(vuzp2_u8(A,A),vuzp1_u8(A,A)),(uint8x16_t)vdupq_n_u64(~2ULL/254)));,auto b=_mm256_set1_epi64x(~2ULL/254);_mm256_storeu_si256(($t*)p,M$(sub,M$(set1,48),M$(cmpeq,_mm256_and_si256(M$(shuffle,_mm256_set1_epi32(x),_mm256_set_epi64x(0,C,C*2,C*3)),b),b)));)p+=$z;}S=p;}$T>V&$O<<(const T&K){F(K);$r}V&$O<<(V&(*A)(V&)){$R A(*this);}};struct W{$T>W&$O<<(const T&K){$r}W&$O<<(W&(*A)(W&)){$R A(*this);}};}namespace std{$f::U i$;$f::V o$;$f::W e$;$f::U&getline($f::U&B,string&K){$f::D A{K};$R B>>A;}$f::V&flush($f::V&B){if(!i$.A.S)B.flush();$R B;}$f::V&endl($f::V&B){$R B<<'\n'<<flush;}$f::W&endl($f::W&B){$R B;}$f::W&flush($f::W&B){$R B;}}$w(LONG WINAPI $x(_EXCEPTION_POINTERS*A){auto C=A->ExceptionRecord;auto B=C->ExceptionInformation[1];if(C->ExceptionCode==2147483649&&B-(ULONG_PTR)std::o$.D<0x40000000){E$(VirtualAlloc((char*)B,16777216,4096,4)&&VirtualAlloc((char*)(B+16777216),4096,4096,260))$R-1;}$R 0;},)
#define freopen(...)if(freopen(__VA_ARGS__)==stdin)std::i$=$f::U{}
#define cin i$
#define cout o$
#ifdef ONLINE_JUDGE
#define cerr e$
#define clog e$
#endif
#line 1 "cp-algo/tree/hld.hpp"
#line 1 "cp-algo/tree/ascending_dfs.hpp"
#line 1 "cp-algo/graph/base.hpp"
#line 1 "cp-algo/graph/edge_types.hpp"
#line 4 "cp-algo/graph/edge_types.hpp"
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
#line 7 "verify/tree/jump.test.cpp"
#include <bits/stdc++.h>
using namespace std;
using namespace cp_algo::graph;
void solve(){int n,q;
cin>>n>>q;
graph g(n);
g.read_edges(n-1);
heavy_light hld(g);
for(int i=0;i<q;i++){node_index u,v;
int k;
cin>>u>>v>>k;
auto ores=hld.jump(u,v,k);
if(!ores){cout<<-1<<'\n';}else{cout<<*ores<<'\n';}}}
signed main(){ios::sync_with_stdio(0);
cin.tie(0);
int t=1;
while(t--){solve();}}