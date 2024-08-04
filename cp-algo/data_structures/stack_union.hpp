#ifndef CP_ALGO_DATA_STRUCTURES_STACK_UNION_HPP
#define CP_ALGO_DATA_STRUCTURES_STACK_UNION_HPP
#include <cstddef>
#include <vector>
namespace cp_algo::data_structures {
    template<class datatype>
    struct stack_union {
        stack_union(int n = 0): head(n), next(1), data(1) {}

        void push(int v, datatype const& vdata) {
            next.push_back(head[v]);
            head[v] = std::size(next) - 1;
            data.push_back(vdata);
        }
        template<typename... Args>
        void emplace(int v, Args&&... vdata) {
            next.push_back(head[v]);
            head[v] = size(next) - 1;
            data.emplace_back(std::forward<Args...>(vdata...));
        }

        size_t size() const {return std::size(head);}
        size_t nodes() const {return std::size(data);}

        std::vector<int> head, next;
        std::vector<datatype> data;
    };
}
#endif // CP_ALGO_DATA_STRUCTURES_STACK_UNION_HPP
