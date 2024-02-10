#ifndef CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP
#define CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP
#include <functional>
#include <algorithm>
#include <cstdint>
namespace cp_algo::data_structures::segment_tree::metas {
    template<typename derived_meta>
    struct base_meta {
        using meta = derived_meta;
        virtual void pull(meta const&, meta const&, int, int) = 0;
        virtual void push(meta*, meta*, int, int) = 0;
    };
}
#endif // CP_ALGO_DATA_STRUCTURES_SEGMENT_TREE_METAS_BASE_HPP