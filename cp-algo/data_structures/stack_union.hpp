#ifndef CP_ALGO_DATA_STRUCTURES_STACK_UNION_HPP
#define CP_ALGO_DATA_STRUCTURES_STACK_UNION_HPP
#include <cstddef>
#include <vector>
namespace cp_algo::data_structures {
    template<class datatype>
    struct stack_union {
        stack_union(size_t n): n(n), head(n), next(1), data(1) {}

        void push(size_t v, datatype const& vdata) {
            next.push_back(head[v]);
            head[v] = size(next) - 1;
            data.push_back(vdata);
        }
        template<typename... Args>
        void emplace(size_t v, Args&&... vdata) {
            next.push_back(head[v]);
            head[v] = size(next) - 1;
            data.emplace_back(std::forward<Args...>(vdata...));
        }

        size_t n;
        std::vector<size_t> head, next;
        std::vector<datatype> data;
    };
}
#endif // CP_ALGO_DATA_STRUCTURES_STACK_UNION_HPP
