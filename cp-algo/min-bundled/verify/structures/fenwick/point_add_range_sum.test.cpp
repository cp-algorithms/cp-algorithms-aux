#line 1 "verify/structures/fenwick/point_add_range_sum.test.cpp"
#define PROBLEM "https:#pragma GCC optimize("Ofast,unroll-loops")
#line 1 "cp-algo/structures/fenwick.hpp"
#include <cassert>
#include <vector>
namespace cp_algo::structures{template<typename Op>
struct inverse_op{};
template<typename T>
struct inverse_op<std::plus<T>>{static T apply(T const&a,T const&b){return a-b;}};
template<typename T>
struct inverse_op<std::multiplies<T>>{static T apply(T const&a,T const&b){return a/b;}};
template<typename T,std::ranges::range Container=std::vector<T>,typename Op=std::plus<T>>
struct fenwick{Op op;
size_t n;
Container data;
fenwick(auto&&range,Op&&op=Op{}):op(std::move(op)){assign(std::move(range));}
void to_prefix_folds(){for(size_t i=1;i<n;i++){if(i+(i&-i)<=n){data[i+(i&-i)]=op(data[i+(i&-i)],data[i]);}}}
void assign(auto&&range){n=size(range)-1;
data=std::move(range);
to_prefix_folds();}
void update(size_t x,T const&v){for(++x;x<=n;x+=x&-x){data[x]=op(data[x],v);}}
T prefix_fold(size_t r)const{assert(r<=n);
T res={};
for(;r;r-=r&-r){res=op(res,data[r]);}
return res;}
T range_fold(size_t l,size_t r)const{return inverse_op<Op>::apply(prefix_fold(r),prefix_fold(l));}
auto prefix_lower_bound(T k)const{size_t x=0;
T pref={};
for(size_t i=std::bit_floor(n);i;i/=2){if(x+i<=n&&op(pref,data[x+i])<=k){pref=op(pref,data[x+i]);
x+=i;}}
return std::pair{x,pref};}};
template<std::ranges::range Container,typename Op>
fenwick(Container&&,Op&&)->fenwick<std::ranges::range_value_t<Container>,Container,Op>;
template<std::ranges::range Container>
fenwick(Container&&)->fenwick<std::ranges::range_value_t<Container>,Container>;
auto maxer=[](auto const&a,auto const&b){return std::max(a,b);};
template<typename T,std::ranges::range Container=std::vector<T>>
struct fenwick_max:fenwick<T,Container,decltype(maxer)>{using fenwick<T,Container,decltype(maxer)>::fenwick;};
template<std::ranges::range Container>
fenwick_max(Container&&)->fenwick_max<std::ranges::range_value_t<Container>,Container>;}
#line 5 "verify/structures/fenwick/point_add_range_sum.test.cpp"
#include <bits/stdc++.h>
using namespace std;
void solve(){int n,q;
cin>>n>>q;
vector<int64_t>a(n+1);
for(auto&it:a|views::drop(1)){cin>>it;}
cp_algo::structures::fenwick me(move(a));
for(int i=0;i<q;i++){int t,x,y;
cin>>t>>x>>y;
if(t==0){me.update(x,y);}else{cout<<me.range_fold(x,y)<<'\n';}}}
signed main(){ios::sync_with_stdio(0);
cin.tie(0);
int t;
t=1;
while(t--){solve();}}