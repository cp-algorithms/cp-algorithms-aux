#ifndef CP_ALGO_STRUCTURES_CSR_HPP
#define CP_ALGO_STRUCTURES_CSR_HPP
#include "../util/big_alloc.hpp"
#include <algorithm>
#include <numeric>
#include <ranges>
#include <span>
namespace cp_algo::structures{template<class datatype>
struct csr{csr():head{0}{}
void new_row(){head.push_back(head.back());}
void push(datatype const&value){data.push_back(value);
head.back()++;}
template<typename... Args>
void emplace(Args&&... args){data.emplace_back(std::forward<Args>(args)...);
head.back()++;}
void reserve_rows(int n){head.reserve(n+1);}
void reserve_data(size_t n){data.reserve(n);}
void reverse_rows(){std::ranges::reverse(data);
std::adjacent_difference(head.begin(),head.end(),head.begin());
std::ranges::reverse(head|std::views::drop(1));
std::partial_sum(head.begin(),head.end(),head.begin());}
size_t size()const{return head.size()-1;}
size_t total_size()const{return data.size();}
auto operator[](this auto&&self,auto row){return std::span(self.data).subspan(self.head[row],self.head[row+1]-self.head[row]);}
auto rows(this auto&&self){return std::views::iota(size_t(0),self.size())
|std::views::transform([&self](auto i){return self[i];});}
big_vector<int>head;
big_vector<datatype>data;};}
#endif