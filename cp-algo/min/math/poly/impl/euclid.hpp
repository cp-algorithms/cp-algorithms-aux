#ifndef CP_ALGO_MATH_POLY_IMPL_EUCLID_HPP
#define CP_ALGO_MATH_POLY_IMPL_EUCLID_HPP
#include "../../affine.hpp"
#include "../../fft.hpp"
#include <functional>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <vector>
#include <tuple>
#include <list>
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math::poly::impl{template<typename poly>using gcd_result=std::pair<std::list<std::decay_t<poly>>,linfrac<std::decay_t<poly>>>;template<typename poly>gcd_result<poly>half_gcd(poly&&A,poly&&B){assert(A.deg()>=B.deg());size_t m=size(A.a)/2;if(B.deg()<(int)m){return{};}auto[ai,R]=A.divmod(B);std::tie(A,B)={B,R};std::list a={ai};auto T=-linfrac(ai).adj();auto advance=[&](size_t k){auto[ak,Tk]=half_gcd(A.div_xk(k),B.div_xk(k));a.splice(end(a),ak);T.prepend(Tk);return Tk;};advance(m).apply(A,B);if constexpr(std::is_reference_v<poly>){advance(2*m-A.deg()).apply(A,B);}else{advance(2*m-A.deg());}return{std::move(a),std::move(T)};}template<typename poly>gcd_result<poly>full_gcd(poly&&A,poly&&B){using poly_t=std::decay_t<poly>;std::list<poly_t>ak;big_vector<linfrac<poly_t>>trs;while(!B.is_zero()){auto[a0,R]=A.divmod(B);ak.push_back(a0);trs.push_back(-linfrac(a0).adj());std::tie(A,B)={B,R};auto[a,Tr]=half_gcd(A,B);ak.splice(end(ak),a);trs.push_back(Tr);}return{ak,std::accumulate(rbegin(trs),rend(trs),linfrac<poly_t>{},std::multiplies{})};}auto convergent(auto L,auto R){using poly=decltype(L)::value_type;if(R==next(L)){return linfrac(*L);}else{int s=std::transform_reduce(L,R,0,std::plus{},std::mem_fn(&poly::deg));auto M=L;for(int c=M->deg();2*c<=s;M++){c+=next(M)->deg();}return convergent(L,M)*convergent(M,R);}}template<typename poly>poly min_rec(poly const&p,size_t d){auto R2=p.mod_xk(d).reversed(d),R1=poly::xk(d);if(R2.is_zero()){return poly(1);}auto[a,Tr]=full_gcd(R1,R2);a.emplace_back();auto pref=begin(a);for(int delta=(int)d-a.front().deg();delta>=0;pref++){delta-=pref->deg()+next(pref)->deg();}return convergent(begin(a),pref).a;}template<typename poly>std::optional<poly>inv_mod(poly p,poly q){assert(!q.is_zero());auto[a,Tr]=full_gcd(q,p);if(q.deg()!=0){return std::nullopt;}return Tr.b/q[0];}}
#pragma GCC pop_options
#endif
